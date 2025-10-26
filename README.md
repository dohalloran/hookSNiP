# hookSNiP 🧬  
_A local web app for allele-specific primer design for hookworm_

<p align="center">
  <a href="#license"><img alt="License: MIT" src="https://img.shields.io/badge/License-MIT-blue.svg"></a>
  <img alt="Python" src="https://img.shields.io/badge/Python-3.x-informational">
  <img alt="Framework" src="https://img.shields.io/badge/FastAPI-server-brightgreen">
  <img alt="PRs welcome" src="https://img.shields.io/badge/PRs-welcome-success">
</p>

> **What is it?** `hookSNiP` runs locally and helps you design allele-specific primers against hookworm genomes.

<p align="center">
  <!-- Replace with a real screenshot or GIF -->
  ![screenshot](https://github.com/dohalloran/hookSNiP/blob/main/image.png?raw=true) 
</p>
---

## ✨ Features
- Local, private processing — bring your own genomes (FASTA).
- Point the app at your genomes with a single env var: `HOOKSNIP_GENOMES_DIR`.
- Lightweight web UI served from your machine.
- No cloud dependencies required.

---

## 🚀 Quickstart

### 1) Create a virtualenv & install
```bash
python3 -m venv .venv && . .venv/bin/activate
pip3 install -r requirements.txt

mkdir -p data/genomes && \
curl -L --fail --retry 3 -o data/genomes/ac_prjna72585.fa.gz "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ancylostoma_caninum/PRJNA72585/ancylostoma_caninum.PRJNA72585.WBPS19.genomic.fa.gz" && \
curl -L --fail --retry 3 -o data/genomes/acey_prjna72583.fa.gz "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ancylostoma_ceylanicum/PRJNA72583/ancylostoma_ceylanicum.PRJNA72583.WBPS19.genomic.fa.gz" && \
curl -L --fail --retry 3 -o data/genomes/acey_prjna231479.fa.gz "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/ancylostoma_ceylanicum/PRJNA231479/ancylostoma_ceylanicum.PRJNA231479.WBPS19.genomic.fa.gz" && \
gunzip -f data/genomes/*.fa.gz

export HOOKSNIP_GENOMES_DIR=/absolute/path/to/genomes

