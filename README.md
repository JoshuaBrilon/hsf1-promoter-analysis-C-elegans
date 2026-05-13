# HSF1 Promoter Analysis in *C. elegans*

Genome-wide analysis of HSF-1 binding motifs in promoter regions of *C. elegans* using a custom computational genomics pipeline.

---

## Overview

This project identifies candidate HSF-1 regulatory targets by scanning promoter regions for Heat Shock Elements (HSEs) and analyzing motif enrichment using FIMO from the MEME Suite.

---

## Repository Structure

```
├── data/
│   ├── raw/               # Input files (not tracked in git)
│   ├── processed/         # Intermediate pipeline outputs
│   └── motif/             # MA2169.1.meme (HSF-1 motif from JASPAR)
├── results/
│   ├── fimo_out_dedupes/  # FIMO output
│   └── plots/             # Histograms and visualizations
├── scripts/               # Pipeline scripts 01–06
└── repo_outputs/          # Key outputs included in this repo
```

---

## Pipeline

Run scripts in order from the repo root. Two steps require external tools between scripts.

**1. Scan genome for HSE motifs**
```
python scripts/01_find_HSE_motifs.py
```

**2. Map motifs to promoter regions**
```
python scripts/02_hsf1_promoter_overlap_analysis.py
```

**3. Annotate promoters with gene names**
```
python scripts/03_annotate_promoters_with_gene_names.py
```

Then run BEDTools:
```bash
bedtools intersect \
  -a data/processed/allHSEs_bed6.txt \
  -b data/processed/promoters_with_gene_names.txt \
  -wa -wb -s > data/processed/promoter_overlaps_with_all_genes.txt
```

**4. Extract promoter sequences**
```
python scripts/04_prepare_promoter_fasta_for_fimo.py
```

**5. Deduplicate promoter regions**
```
python scripts/05_deduplicate_fimo_promoter_hits.py
```

Then build the background model and run FIMO:
```bash
fasta-get-markov -m 1 \
  data/processed/allPromoters_dedup.fasta \
  data/processed/background_dedup.txt

fimo --oc results/fimo_out_dedupes \
  --bgfile data/processed/background_dedup.txt \
  data/motif/MA2169.1.meme \
  data/processed/allPromoters_dedup.fasta
```

**6. Summarize FIMO results by gene**
```
python scripts/06_summarize_fimo_by_gene.py
```

---

## Outputs

- **fimo_gene_summary.txt** — genes ranked by HSF-1 motif enrichment
- **promoters_with_gene_names_dedup.txt** — deduplicated promoter dataset
- **allPromoters_dedup.fasta** — promoter sequences used as FIMO input

---

## Key Result

Deduplicating promoter regions reduced motif counts by ~19%, indicating that redundancy from multi-isoform gene models can inflate motif-based predictions.

---

## Data Sources

| File | Description | Source |
|------|-------------|--------|
| `libuda_N2_genome.fasta` | *C. elegans* N2 reference genome | Libuda Lab (not publicly distributed) |
| `N2.genome.annotations.gff3` | Genome annotations | Libuda Lab |
| `C_elegans.current.geneIDs.txt` | WormBase gene ID table | [WormBase](https://wormbase.org) |
| `MA2169.1.meme` | HSF-1 binding motif (JASPAR ID: MA2169.1) | [JASPAR 2024](https://jaspar.elixir.no/matrix/MA2169.1/) |

> Raw genome and annotation files are not included in this repository.  
> Contact the Libuda Lab at the University of Oregon for access.

---

## Requirements

**Python packages** (Python ≥ 3.10)

```
pip install biopython pandas matplotlib seaborn
```

**External tools**
- [MEME Suite](https://meme-suite.org) — for `fimo` and `fasta-get-markov`
- [BEDTools](https://bedtools.readthedocs.io) — for promoter overlap

---

## Author

Joshua Brilon  
University of Oregon
