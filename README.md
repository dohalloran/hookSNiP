
# hookSNiP 

hookSNiP: allele-specific primer design for hookworm

## Quickstart

```bash
python3 -m venv .venv && . .venv/bin/activate
pip3 install -r requirements.txt
mkdir -p data/genomes && \
curl -L --fail --retry 3 -o data/genomes/ac_prjna72585.fa.gz "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ancylostoma_caninum/PRJNA72585/ancylostoma_caninum.PRJNA72585.WBPS19.genomic.fa.gz" && \
curl -L --fail --retry 3 -o data/genomes/acey_prjna72583.fa.gz "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ancylostoma_ceylanicum/PRJNA72583/ancylostoma_ceylanicum.PRJNA72583.WBPS19.genomic.fa.gz" && \
curl -L --fail --retry 3 -o data/genomes/acey_prjna231479.fa.gz "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ancylostoma_ceylanicum/PRJNA231479/ancylostoma_ceylanicum.PRJNA231479.WBPS19.genomic.fa.gz" && \
gunzip -f data/genomes/*.fa.gz
# copy your built-in genomes to data/genomes (plain .fa or .fa.gz)
uvicorn server.main:app --reload --port 8000
# open http://localhost:8000/
```

## Built-in genomes
Place files under `data/genomes` or set `HOOKSNIP_GENOMES_DIR=/absolute/path/to/genomes`.

Accepted:
- `.fa`, `.fasta`, `.fna` (plain text)
- `.fa.gz` (server will decompress into a temp workspace)


