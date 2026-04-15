import os
import os.path as osp
import sys
import json
import random
import html as html_mod

# Where the server tells us to write (set per request in main.py)
OUTDIR = os.environ.get("OUTDIR", "./out")

# ---------------------------------------------------------------------------
# FIX #1: The original used a builtins.open monkey-patch to redirect "/out/…"
# paths written by run_pipeline to OUTDIR. That monkey-patch is fragile
# (affects all open() calls in the process, including stdlib internals), and
# it was set at module-import time so OUTDIR was frozen to the env-var value
# rather than the per-request value set by main.py.
#
# The fix: run_pipeline now accepts outdir as an explicit parameter, and all
# output file writes use that parameter directly. No monkey-patching needed.
# ---------------------------------------------------------------------------


def read_fasta(path):
    seqs = {}
    name = None
    parts = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = "".join(parts).upper()
    return seqs


def parse_vcf_snvs(path):
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            chrom, pos, vid, ref, alts = p[:5]
            try:
                pos = int(pos)
            except ValueError:
                continue
            if len(ref) != 1:
                continue
            alt = alts.split(",")[0]
            if len(alt) != 1:
                continue
            rows.append(
                dict(chrom=chrom, pos=pos, id=vid, ref=ref.upper(), alt=alt.upper())
            )
    return rows


def parse_gff_strand_index(path):
    keep = {"gene", "mRNA", "transcript", "CDS", "exon"}
    idx = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if not ln or ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            typ = p[2]
            if typ not in keep:
                continue
            try:
                s = int(p[3])
                e = int(p[4])
            except ValueError:
                continue
            st = p[6]
            if st not in ("+", "-"):
                continue
            idx.setdefault(p[0], []).append((s, e, st))
    for c in idx:
        idx[c].sort()
    return idx


def strand_for(idx, chrom, pos):
    arr = idx.get(chrom)
    if not arr:
        return "+"
    for s, e, st in arr:
        if s <= pos <= e:
            return st
    return "+"


_comp = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(_comp)[::-1]


def tm(s):
    s = s.upper()
    return 2 * (s.count("A") + s.count("T")) + 4 * (s.count("G") + s.count("C"))


def gc(s):
    n = max(1, len(s))
    s = s.upper()
    return 100.0 * (s.count("G") + s.count("C")) / n


def hpoly(s):
    if not s:
        return 0
    best = cur = 1
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


DEFAULTS = dict(
    len_lo=18,
    len_hi=25,
    tm_target=60.0,
    tm_tol=4.0,
    gc_lo=35.0,
    gc_hi=65.0,
    amp_lo=120,
    amp_hi=250,
    homopoly_max=5,
)


def fetch_window(seqs, chrom, pos, flank):
    seq = seqs.get(chrom)
    if seq is None:
        raise RuntimeError(f"Contig not in FASTA: {chrom}")
    i = pos - 1
    if i < 0 or i >= len(seq):
        raise RuntimeError(f"POS out of range on {chrom}")
    up = seq[max(0, i - flank) : i]
    dn = seq[i + 1 : min(len(seq), i + 1 + flank)]
    return up, dn, seq[i]


def mm_base(b, avoid):
    for x in ("A", "T"):
        if x != b and x != avoid:
            return x
    for x in ("A", "T"):
        if x != b:
            return x
    return "A"


def fmt_left(seq, mm_applied):
    # FIX #2: Primer sequences must be HTML-escaped before insertion into HTML
    # to avoid XSS if a sequence ever contains '<', '>', or '&' (unlikely for
    # DNA but still a correctness issue — and the span tags must survive escaping).
    if len(seq) >= 3:
        body = [html_mod.escape(c) for c in seq]
        if mm_applied:
            body[-3] = f"<span class='mm'>{body[-3]}</span>"
        body[-1] = f"<span class='alt'>{body[-1]}</span>"
        return "<code>" + "".join(body) + "</code>"
    elif len(seq) >= 1:
        escaped = html_mod.escape(seq)
        return "<code>" + escaped[:-1] + "<span class='alt'>" + escaped[-1] + "</span></code>"
    else:
        return "<code></code>"


def score_pair(Lseq, Rseq):
    d = abs(tm(Lseq) - DEFAULTS["tm_target"]) + abs(tm(Rseq) - DEFAULTS["tm_target"])
    gc_pen = 0.0
    for s in (Lseq, Rseq):
        g = gc(s)
        if g < DEFAULTS["gc_lo"] or g > DEFAULTS["gc_hi"]:
            gc_pen += 10.0
    return max(0.0, 100.0 - (d + gc_pen))


def design_one(seqs, gff_idx, rec, flank, top_per_allele, insert_mm):
    chrom, pos, vid, ref, alt = (
        rec["chrom"],
        rec["pos"],
        rec["id"],
        rec["ref"],
        rec["alt"],
    )
    up, dn, ref_base = fetch_window(seqs, chrom, pos, flank)
    if ref_base.upper() != ref.upper():
        return dict(chrom=chrom, pos=pos, id=vid, strand="?", alleles=[], note="REF mismatch")
    st = strand_for(gff_idx, chrom, pos)
    out = []
    for label, allele in (("REF", ref), ("ALT", alt)):
        Lc = []
        Rc = []
        # LEFT (allele-specific)
        for L in range(DEFAULTS["len_lo"], DEFAULTS["len_hi"] + 1):
            if st == "+":
                if len(up) < L - 1:
                    continue
                left = up[-(L - 1) :] + allele
            else:
                if len(dn) < L - 1:
                    continue
                proto = allele + dn[: L - 1]
                left = rc(proto)
            mm_applied = False
            if (label == "ALT") and insert_mm and L >= 3:
                idx_mm = -3
                nb = mm_base(left[idx_mm], left[-1])
                left = left[:idx_mm] + nb + left[idx_mm + 1 :]
                mm_applied = True
            if hpoly(left) > DEFAULTS["homopoly_max"]:
                continue
            if abs(tm(left) - DEFAULTS["tm_target"]) > DEFAULTS["tm_tol"]:
                continue
            if not (DEFAULTS["gc_lo"] <= gc(left) <= DEFAULTS["gc_hi"]):
                continue
            Lc.append((left, L, mm_applied))
        # RIGHT (non-discriminating)
        if st == "+":
            for LR in range(DEFAULTS["len_lo"], DEFAULTS["len_hi"] + 1):
                if len(dn) < LR:
                    continue
                for k in range(0, len(dn) - LR + 1):
                    r = rc(dn[k : k + LR])
                    if hpoly(r) > DEFAULTS["homopoly_max"]:
                        continue
                    if abs(tm(r) - DEFAULTS["tm_target"]) > DEFAULTS["tm_tol"]:
                        continue
                    if not (DEFAULTS["gc_lo"] <= gc(r) <= DEFAULTS["gc_hi"]):
                        continue
                    Rc.append((k, r, LR))
        else:
            for LR in range(DEFAULTS["len_lo"], DEFAULTS["len_hi"] + 1):
                if len(up) < LR:
                    continue
                for k in range(0, len(up) - LR + 1):
                    end = len(up) - k
                    start = end - LR
                    if start < 0:
                        continue
                    r = up[start:end]
                    if hpoly(r) > DEFAULTS["homopoly_max"]:
                        continue
                    if abs(tm(r) - DEFAULTS["tm_target"]) > DEFAULTS["tm_tol"]:
                        continue
                    if not (DEFAULTS["gc_lo"] <= gc(r) <= DEFAULTS["gc_hi"]):
                        continue
                    Rc.append((k, r, LR))
        pairs = []
        for Lseq, LL, mm_ap in Lc:
            for k, Rseq, LR in Rc:
                amplen = LL + k + LR
                if not (DEFAULTS["amp_lo"] <= amplen <= DEFAULTS["amp_hi"]):
                    continue
                pairs.append(
                    dict(
                        left=Lseq,
                        right=Rseq,
                        left_len=LL,
                        right_len=LR,
                        left_tm=tm(Lseq),
                        right_tm=tm(Rseq),
                        left_gc=gc(Lseq),
                        right_gc=gc(Rseq),
                        amplicon=amplen,
                        score=score_pair(Lseq, Rseq),
                        left_fmt=fmt_left(Lseq, mm_ap),
                    )
                )
        pairs.sort(key=lambda p: (-p["score"], abs(p["amplicon"] - 180)))
        out.append(dict(allele=label, pairs=pairs[: max(1, int(top_per_allele))]))
    return dict(chrom=chrom, pos=pos, id=vid, strand=st, alleles=out)


def quick_check(seqs, rows):
    contigs = set(seqs.keys())
    seen = set(r["chrom"] for r in rows)
    missing = sorted(c for c in seen if c not in contigs)
    if missing:
        # FIX #3: Cap how many missing contigs are shown so the error message
        # doesn't become enormous with large VCFs.
        shown = missing[:10]
        suffix = f" … and {len(missing)-10} more" if len(missing) > 10 else ""
        raise RuntimeError(
            "Contigs in VCF not present in FASTA: " + ", ".join(shown) + suffix
        )
    # FIX #4 (carried from previous session): Removed the hard REF-mismatch
    # check that aborted the entire run if any sampled variant mismatched.
    # design_one already handles per-variant mismatches gracefully with a note.


def render_table(results):
    P = []
    P.append("<h2>Results</h2>")
    P.append(
        "<table><thead><tr>"
        "<th>Variant</th><th>Allele</th><th>Strand</th>"
        "<th>LEFT (5\u2032\u21923\u2032)</th><th>RIGHT (5\u2032\u21923\u2032)</th>"
        "<th>Len (L/R)</th><th>Tm (L/R)</th><th>GC% (L/R)</th>"
        "<th>Amplicon</th><th>Score</th>"
        "</tr></thead><tbody>"
    )
    for R in results:
        note = R.get("note", "")
        note_html = f" \u2014 {html_mod.escape(note)}" if note else ""
        P.append(
            f"<tr><th colspan=10 style='text-align:left;background:#f6f8fa'>"
            f"{html_mod.escape(R['chrom'])}:{R['pos']} (strand {R['strand']}){note_html}"
            f"</th></tr>"
        )
        for A in R["alleles"]:
            if not A["pairs"]:
                P.append(
                    f"<tr><td>{html_mod.escape(R['chrom'])}:{R['pos']}</td>"
                    f"<td>{A['allele']}</td><td>{R['strand']}</td>"
                    f"<td colspan=7 style='color:#777'>No pairs under constraints</td></tr>"
                )
                continue
            for p in A["pairs"]:
                # FIX #5: right primer sequence was not HTML-escaped in render_table
                right_escaped = html_mod.escape(p["right"])
                P.append(
                    "<tr>"
                    f"<td>{html_mod.escape(R['chrom'])}:{R['pos']}</td>"
                    f"<td>{A['allele']}</td>"
                    f"<td>{R['strand']}</td>"
                    f"<td>{p['left_fmt']}</td>"
                    f"<td><code>{right_escaped}</code></td>"
                    f"<td>{p['left_len']}/{p['right_len']}</td>"
                    f"<td>{int(p['left_tm'])}/{int(p['right_tm'])}</td>"
                    f"<td>{p['left_gc']:.1f}/{p['right_gc']:.1f}</td>"
                    f"<td>{p['amplicon']}</td>"
                    f"<td>{p['score']:.1f}</td>"
                    "</tr>"
                )
    P.append("</tbody></table>")
    P.append(
        "<div class='legend'>Legend: "
        "<span class='mm'>blue</span>=\u22122 mismatch (LEFT \u22122); "
        "<span class='alt'>red</span>=3\u2032 allele base.</div>"
    )
    return "".join(P)


# FIX #1 (cont): outdir is now an explicit parameter, not a global/env-var.
def run_pipeline(vcf_path, fa_path, gff_path, flank, top_per_allele, insert_mismatch, outdir=None):
    if outdir is None:
        outdir = OUTDIR
    os.makedirs(outdir, exist_ok=True)

    seqs = read_fasta(fa_path)
    rows = parse_vcf_snvs(vcf_path)
    gff = parse_gff_strand_index(gff_path)
    quick_check(seqs, rows)
    results = []
    processed = 0
    for rec in rows:
        results.append(design_one(seqs, gff, rec, flank, top_per_allele, insert_mismatch))
        processed += 1

    head = (
        "<h2>Run</h2><table>"
        f"<tr><th>VCF</th><td>{html_mod.escape(osp.basename(vcf_path))}</td></tr>"
        f"<tr><th>FASTA</th><td>{html_mod.escape(osp.basename(fa_path))}</td></tr>"
        f"<tr><th>GFF3</th><td>{html_mod.escape(osp.basename(gff_path))}</td></tr>"
        f"<tr><th>Flank</th><td>{flank}</td></tr>"
        f"<tr><th>Top per allele</th><td>{top_per_allele}</td></tr>"
        f"<tr><th>\u22122 mismatch</th><td>{'enabled' if insert_mismatch else 'disabled'}</td></tr>"
        f"<tr><th>Variants processed</th><td>{processed}</td></tr>"
        "</table>"
    )
    tbl = render_table(results)

    primers_path = osp.join(outdir, "primers.html")
    table_path = osp.join(outdir, "table.html")

    with open(primers_path, "w", encoding="utf-8") as fh:
        fh.write("<!doctype html><meta charset='utf-8'><title>primers</title>" + head + tbl)
    with open(table_path, "w", encoding="utf-8") as fh:
        fh.write(tbl)

    return json.dumps(dict(ok=True, processed=processed))
