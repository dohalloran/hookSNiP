import os, builtins, os.path as osp

# Where the server tells us to write (set per request in main.py)
OUTDIR = os.environ.get("OUTDIR", "./out")

# --- Redirect any use of "/out" to OUTDIR ---
_real_open = builtins.open
def _open_proxy(file, *args, **kwargs):
    if isinstance(file, (str, os.PathLike)):
        p = os.fspath(file)
        if p == "/out" or p.startswith("/out/"):
            # map "/out/foo" -> f"{OUTDIR}/foo"
            mapped = osp.join(OUTDIR, p[5:] if p != "/out" else "")
            os.makedirs(osp.dirname(mapped) or OUTDIR, exist_ok=True)
            file = mapped
    return _real_open(file, *args, **kwargs)
builtins.open = _open_proxy

# Also catch explicit makedirs("/out/...") calls
_real_makedirs = os.makedirs
def _makedirs_proxy(path, *args, **kwargs):
    if isinstance(path, (str, os.PathLike)):
        p = os.fspath(path)
        if p == "/out" or p.startswith("/out/"):
            path = osp.join(OUTDIR, p[5:] if p != "/out" else "")
    return _real_makedirs(path, *args, **kwargs)
os.makedirs = _makedirs_proxy




import os
OUTDIR = os.environ.get("OUTDIR", "./out")

import sys, os, json, random

def read_fasta(path):
    seqs = {}; name=None; parts=[]
    with open(path,'r',encoding='utf-8',errors='ignore') as fh:
        for line in fh:
            if not line: continue
            if line.startswith('>'):
                if name is not None: seqs[name] = ''.join(parts).upper()
                name = line[1:].strip().split()[0]; parts=[]
            else:
                parts.append(line.strip())
    if name is not None: seqs[name] = ''.join(parts).upper()
    return seqs

def parse_vcf_snvs(path):
    rows=[]
    with open(path,'r',encoding='utf-8',errors='ignore') as fh:
        for ln in fh:
            if not ln or ln.startswith('#'): continue
            p = ln.rstrip('\n').split('\t')
            if len(p) < 5: continue
            chrom, pos, vid, ref, alts = p[:5]
            try: pos = int(pos)
            except: continue
            if len(ref) != 1: continue
            alt = alts.split(',')[0]
            if len(alt) != 1: continue
            rows.append(dict(chrom=chrom,pos=pos,id=vid,ref=ref.upper(),alt=alt.upper()))
    return rows

def parse_gff_strand_index(path):
    keep = {'gene','mRNA','transcript','CDS','exon'}
    idx={}
    with open(path,'r',encoding='utf-8',errors='ignore') as fh:
        for ln in fh:
            if not ln or ln.startswith('#'): continue
            p=ln.rstrip('\n').split('\t')
            if len(p) < 9: continue
            typ=p[2]
            if typ not in keep: continue
            try:
                s=int(p[3]); e=int(p[4])
            except: continue
            st=p[6]
            if st not in ('+','-'): continue
            idx.setdefault(p[0], []).append((s,e,st))
    for c in idx: idx[c].sort()
    return idx

def strand_for(idx, chrom, pos):
    arr = idx.get(chrom)
    if not arr: return '+'
    for s,e,st in arr:
        if s <= pos <= e: return st
    return '+'

_comp = str.maketrans('ACGTN','TGCAN')
def rc(s): return s.translate(_comp)[::-1]
def tm(s): s=s.upper(); return 2*(s.count('A')+s.count('T')) + 4*(s.count('G')+s.count('C'))
def gc(s): n=max(1,len(s)); s=s.upper(); return 100.0*(s.count('G')+s.count('C'))/n
def hpoly(s):
    if not s: return 0
    best=cur=1
    for i in range(1,len(s)):
        if s[i]==s[i-1]: cur+=1; best=max(best,cur)
        else: cur=1
    return best

DEFAULTS = dict(len_lo=18,len_hi=25, tm_target=60.0, tm_tol=4.0, gc_lo=35.0, gc_hi=65.0,
                amp_lo=120, amp_hi=250, homopoly_max=5)

def fetch_window(seqs, chrom, pos, flank):
    seq = seqs.get(chrom)
    if seq is None: raise RuntimeError(f"Contig not in FASTA: {chrom}")
    i=pos-1
    if i<0 or i>=len(seq): raise RuntimeError(f"POS out of range on {chrom}")
    up = seq[max(0,i-flank):i]
    dn = seq[i+1:min(len(seq), i+1+flank)]
    return up, dn, seq[i]

def mm_base(b, avoid):
    # Best: A/T that differ from both the original base and the 3′ base
    for x in ('A','T'):
        if x != b and x != avoid:
            return x
    # Otherwise: A/T that differ from the original base (even if == 3′ base)
    for x in ('A','T'):
        if x != b:
            return x
    # Fallback (edge cases): A
    return 'A'

def fmt_left(seq, mm_applied):
    if len(seq) >= 3:
        body = list(seq)
        if mm_applied: body[-3] = f"<span class='mm'>{body[-3]}</span>"
        body[-1] = f"<span class='alt'>{body[-1]}</span>"
        return "<code>" + "".join(body) + "</code>"
    elif len(seq) >= 1:
        return "<code>"+seq[:-1]+"<span class='alt'>"+seq[-1]+"</span></code>"
    else:
        return "<code>"+seq+"</code>"

def score_pair(Lseq, Rseq):
    d = abs(tm(Lseq)-DEFAULTS['tm_target']) + abs(tm(Rseq)-DEFAULTS['tm_target'])
    gc_pen = 0.0
    for s in (Lseq,Rseq):
        g = gc(s)
        if g < DEFAULTS['gc_lo'] or g > DEFAULTS['gc_hi']: gc_pen += 10.0
    return max(0.0, 100.0 - (d + gc_pen))

def design_one(seqs, gff_idx, rec, flank, top_per_allele, insert_mm):
    chrom,pos,vid,ref,alt = rec['chrom'], rec['pos'], rec['id'], rec['ref'], rec['alt']
    up, dn, ref_base = fetch_window(seqs, chrom, pos, flank)
    if ref_base.upper() != ref.upper():
        return dict(chrom=chrom,pos=pos,id=vid,strand='?',alleles=[],note='REF mismatch')
    st = strand_for(gff_idx, chrom, pos)
    out=[]
    for label, allele in (('REF',ref),('ALT',alt)):
        Lc=[]; Rc=[]
        # LEFT (allele-specific)
        for L in range(DEFAULTS['len_lo'], DEFAULTS['len_hi']+1):
            if st == '+':
                if len(up) < L-1: continue
                left = up[-(L-1):] + allele
            else:
                if len(dn) < L-1: continue
                proto = allele + dn[:L-1]
                left = rc(proto)
            mm_applied=False
            if (label == 'ALT') and insert_mm and L>=3:
                i=-3; nb = mm_base(left[i], left[-1]); left = left[:i] + nb + left[i+1:]; mm_applied=True
            if hpoly(left) > DEFAULTS['homopoly_max']: continue
            if abs(tm(left)-DEFAULTS['tm_target']) > DEFAULTS['tm_tol']: continue
            if not (DEFAULTS['gc_lo'] <= gc(left) <= DEFAULTS['gc_hi']): continue
            Lc.append((left,L,mm_applied))
        # RIGHT (non-discriminating)
        if st == '+':
            # from downstream, reverse-complement segment
            for LR in range(DEFAULTS['len_lo'], DEFAULTS['len_hi']+1):
                if len(dn) < LR: continue
                for k in range(0, len(dn)-LR+1):
                    r = rc(dn[k:k+LR])
                    if hpoly(r) > DEFAULTS['homopoly_max']: continue
                    if abs(tm(r)-DEFAULTS['tm_target']) > DEFAULTS['tm_tol']: continue
                    if not (DEFAULTS['gc_lo'] <= gc(r) <= DEFAULTS['gc_hi']): continue
                    Rc.append((k,r,LR))
        else:
            # from upstream, forward segment
            for LR in range(DEFAULTS['len_lo'], DEFAULTS['len_hi']+1):
                if len(up) < LR: continue
                for k in range(0, len(up)-LR+1):
                    end = len(up)-k
                    start = end-LR
                    if start < 0: continue
                    r = up[start:end]
                    if hpoly(r) > DEFAULTS['homopoly_max']: continue
                    if abs(tm(r)-DEFAULTS['tm_target']) > DEFAULTS['tm_tol']: continue
                    if not (DEFAULTS['gc_lo'] <= gc(r) <= DEFAULTS['gc_hi']): continue
                    Rc.append((k,r,LR))
        pairs=[]
        for (Lseq,LL,mm_ap) in Lc:
            for (k,Rseq,LR) in Rc:
                amplen = LL + k + LR
                if not (DEFAULTS['amp_lo'] <= amplen <= DEFAULTS['amp_hi']): continue
                pairs.append(dict(
                    left=Lseq,right=Rseq,left_len=LL,right_len=LR,
                    left_tm=tm(Lseq),right_tm=tm(Rseq),
                    left_gc=gc(Lseq),right_gc=gc(Rseq),
                    amplicon=amplen,score=score_pair(Lseq,Rseq),
                    left_fmt=fmt_left(Lseq, mm_ap)))
        pairs.sort(key=lambda p:(-p['score'], abs(p['amplicon']-180)))
        out.append(dict(allele=label, pairs=pairs[:max(1,int(top_per_allele))]))
    return dict(chrom=chrom,pos=pos,id=vid,strand=st,alleles=out)

def quick_check(seqs, rows):
    contigs=set(seqs.keys()); seen=set(r['chrom'] for r in rows)
    missing=sorted(c for c in seen if c not in contigs)
    if missing: raise RuntimeError("Contigs in VCF not present in FASTA: " + ", ".join(missing))
    # light REF check
    sample = rows[:]; random.shuffle(sample); sample = sample[:min(50, len(sample))]
    mism=0
    for r in sample:
        s = seqs.get(r['chrom']); i=r['pos']-1
        if s is None or i<0 or i>=len(s) or s[i].upper()!=r['ref'].upper(): mism+=1
    if mism: raise RuntimeError(f"Reference base mismatches for {mism}/{len(sample)} tested sites")

def render_table(results):
    P=[]
    P.append("<h2>Results</h2>")
    P.append("<table><thead><tr><th>Variant</th><th>Allele</th><th>Strand</th><th>LEFT (5′→3′)</th><th>RIGHT (5′→3′)</th><th>Len (L/R)</th><th>Tm (L/R)</th><th>GC% (L/R)</th><th>Amplicon</th><th>Score</th></tr></thead><tbody>")
    for R in results:
        P.append(f"<tr><th colspan=10 style='text-align:left;background:#f6f8fa'>{R['chrom']}:{R['pos']} (strand {R['strand']})</th></tr>")
        for A in R['alleles']:
            if not A['pairs']:
                P.append(f"<tr><td>{R['chrom']}:{R['pos']}</td><td>{A['allele']}</td><td>{R['strand']}</td><td colspan=7 style='color:#777'>No pairs under constraints</td></tr>")
                continue
            for p in A['pairs']:
                P.append("<tr>"
                         f"<td>{R['chrom']}:{R['pos']}</td>"
                         f"<td>{A['allele']}</td>"
                         f"<td>{R['strand']}</td>"
                         f"<td>{p['left_fmt']}</td>"
                         f"<td><code>{p['right']}</code></td>"
                         f"<td>{p['left_len']}/{p['right_len']}</td>"
                         f"<td>{int(p['left_tm'])}/{int(p['right_tm'])}</td>"
                         f"<td>{p['left_gc']:.1f}/{p['right_gc']:.1f}</td>"
                         f"<td>{p['amplicon']}</td>"
                         f"<td>{p['score']:.1f}</td>"
                         "</tr>")
    P.append("</tbody></table>")
    P.append("<div class='legend'>Legend: <span class='mm'>blue</span>=−2 mismatch (LEFT −2); <span class='alt'>red</span>=3′ allele base.</div>")
    return "".join(P)

def run_pipeline(vcf_path, fa_path, gff_path, flank, top_per_allele, insert_mismatch):
    seqs = read_fasta(fa_path)
    rows = parse_vcf_snvs(vcf_path)
    gff = parse_gff_strand_index(gff_path)
    quick_check(seqs, rows)
    results=[]; processed=0
    for rec in rows:
        results.append(design_one(seqs, gff, rec, flank, top_per_allele, insert_mismatch))
        processed += 1
    os.makedirs("/out", exist_ok=True)
    head = ("<h2>Run</h2><table>"
            f"<tr><th>VCF</th><td>{os.path.basename(vcf_path)}</td></tr>"
            f"<tr><th>FASTA</th><td>{os.path.basename(fa_path)}</td></tr>"
            f"<tr><th>GFF3</th><td>{os.path.basename(gff_path)}</td></tr>"
            f"<tr><th>Flank</th><td>{flank}</td></tr>"
            f"<tr><th>Top per allele</th><td>{top_per_allele}</td></tr>"
            f"<tr><th>−2 mismatch</th><td>{'enabled' if insert_mismatch else 'disabled'}</td></tr>"
            f"<tr><th>Variants processed</th><td>{processed}</td></tr></table>")
    tbl = render_table(results)
    with open(f"{OUTDIR}/primers.html","w",encoding="utf-8") as fh:
        fh.write("<!doctype html><meta charset='utf-8'><title>primers</title>"+head+tbl)
    with open(f"{OUTDIR}/table.html","w",encoding="utf-8") as fh:
        fh.write(tbl)
    return json.dumps(dict(ok=True, processed=processed))
