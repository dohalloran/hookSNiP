import os
import tempfile
import shutil
import gzip
import json
from pathlib import Path

from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from . import py_core  # run_pipeline is defined there

app = FastAPI(title="HookSnip Server")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

SERVER_GENOMES_DIR = Path(
    os.environ.get(
        "HOOKSNIP_GENOMES_DIR",
        str(Path(__file__).resolve().parent.parent / "data" / "genomes"),
    )
)


def _save_upload(dst: Path, up: UploadFile):
    with dst.open("wb") as f:
        while True:
            chunk = up.file.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def _ensure_plain_fasta(src_path: Path, dst_path: Path):
    if src_path.suffix.lower() == ".gz":
        with gzip.open(src_path, "rb") as gz, dst_path.open("wb") as out:
            shutil.copyfileobj(gz, out)
    else:
        shutil.copyfile(src_path, dst_path)


# FIX #6: Validate uploaded file extensions server-side. The original accepted
# any file for gff3 and vcf — a user (or bad actor) could upload arbitrary
# content. We now reject clearly wrong file types before running the pipeline.
_FASTA_EXTS = {".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"}
_GFF_EXTS   = {".gff", ".gff3"}
_VCF_EXTS   = {".vcf"}

def _check_ext(filename: str, allowed: set, label: str):
    name_lower = filename.lower()
    if not any(name_lower.endswith(ext) for ext in allowed):
        raise HTTPException(
            status_code=400,
            detail=f"{label} must be one of: {', '.join(sorted(allowed))}. Got: {filename!r}",
        )


@app.post("/api/run")
async def run(
    srcmode: str = Form(...),
    genome_name: str = Form(None),
    flank: int = Form(300),
    top_per_allele: int = Form(1),
    insert_mismatch: bool = Form(False),
    genome: UploadFile = File(None),
    gff3: UploadFile = File(...),
    vcf: UploadFile = File(...),
):
    # FIX #7: Clamp parameters to safe ranges server-side so malformed requests
    # (or a tweaked frontend) can't trigger degenerate pipeline behaviour.
    flank         = max(50,  min(flank,         2000))
    top_per_allele = max(1,  min(top_per_allele, 10))

    # Validate file extensions before doing any disk I/O.
    _check_ext(gff3.filename, _GFF_EXTS, "GFF3 file")
    _check_ext(vcf.filename,  _VCF_EXTS, "VCF file")
    if srcmode == "custom" and genome is not None:
        _check_ext(genome.filename, _FASTA_EXTS, "Custom genome file")

    work = Path(tempfile.mkdtemp(prefix="hooksnip_"))
    out_dir = work / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        genome_path = work / "genome.fa"
        gff_path    = work / "anno.gff3"
        vcf_path    = work / "variants.vcf"

        _save_upload(gff_path, gff3)
        _save_upload(vcf_path, vcf)

        if srcmode == "built-in":
            if not genome_name:
                raise HTTPException(status_code=400, detail="Missing genome_name for built-in mode.")
            # FIX #8: Prevent path traversal — genome_name must be a plain
            # filename with no directory components.
            if "/" in genome_name or "\\" in genome_name or genome_name.startswith("."):
                raise HTTPException(status_code=400, detail="Invalid genome_name.")
            cand_plain = SERVER_GENOMES_DIR / genome_name
            gz_name    = genome_name if genome_name.endswith(".gz") else genome_name + ".gz"
            cand_gz    = SERVER_GENOMES_DIR / gz_name
            if cand_plain.exists():
                _ensure_plain_fasta(cand_plain, genome_path)
            elif cand_gz.exists():
                _ensure_plain_fasta(cand_gz, genome_path)
            else:
                raise HTTPException(
                    status_code=404,
                    detail=f"Built-in genome not found on server: {genome_name}",
                )
        else:
            if genome is None:
                raise HTTPException(status_code=400, detail="Please upload a genome file for custom mode.")
            up_path = work / Path(genome.filename).name   # FIX #9: strip any path prefix from filename
            _save_upload(up_path, genome)
            _ensure_plain_fasta(up_path, genome_path)

        # FIX #1 (cont): pass out_dir explicitly; py_core no longer uses a global.
        info_raw = py_core.run_pipeline(
            str(vcf_path),
            str(genome_path),
            str(gff_path),
            flank=flank,
            top_per_allele=top_per_allele,
            insert_mismatch=insert_mismatch,
            outdir=str(out_dir),
        )
        try:
            info = json.loads(info_raw) if isinstance(info_raw, str) else info_raw
        except Exception:
            info = {"ok": False, "error": "Pipeline returned non-JSON response", "raw": str(info_raw)}

        table_html   = (out_dir / "table.html").read_text(encoding="utf-8")
        primers_html = (out_dir / "primers.html").read_text(encoding="utf-8")

        return JSONResponse(
            {"ok": True, "info": info, "table_html": table_html, "primers_html": primers_html}
        )

    except HTTPException:
        raise
    except Exception as e:
        # FIX #10: Don't leak raw tracebacks/internal paths to the client.
        # Log the full error server-side, return only the message.
        import traceback
        traceback.print_exc()
        return JSONResponse({"ok": False, "error": str(e)}, status_code=500)
    finally:
        shutil.rmtree(work, ignore_errors=True)


PUBLIC_DIR = Path(__file__).resolve().parent.parent / "public"
app.mount("/", StaticFiles(directory=str(PUBLIC_DIR), html=True), name="static")
