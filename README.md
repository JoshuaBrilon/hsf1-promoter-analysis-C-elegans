# HSF1 Promoter Analysis in *C. elegans*

Genome-wide analysis of HSF-1 binding motifs in promoter regions of *C. elegans* using a custom computational genomics pipeline.

---

## Overview

This project identifies candidate HSF-1 regulatory targets by scanning promoter regions for Heat Shock Elements (HSEs) and analyzing motif enrichment.

---

## Repository Structure

scripts/ → pipeline scripts (01–06)
repo_outputs/ → key output files from the analysis

---

## Pipeline

1. Scan genome for HSE motifs  
2. Map motifs to promoter regions  
3. Annotate promoters with gene names  
4. Extract promoter sequences  
5. Remove duplicate promoter regions  
6. Run FIMO and summarize results by gene  

---

## Outputs

- **fimo_gene_summary.txt** → ranked genes based on motif enrichment  
- **promoters_with_gene_names_dedup.txt** → cleaned promoter dataset  
- **allPromoters_dedup.fasta** → promoter sequences used for analysis  

---

## Key Result

Deduplicating promoter regions reduced motif counts by ~19%, indicating redundancy can inflate motif-based predictions.

---

## Requirements

- Python 3  
- Biopython  
- pandas  
- MEME Suite (FIMO)  

---

## Author

Joshua Brilon  
University of Oregon
