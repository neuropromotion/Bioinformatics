# Loss of tissue specificity and recurrent pan-cancer activation define a conserved oncogenic microRNA program

Overview

MicroRNAs are post-transcriptional regulators capable of coordinately modulating large gene networks and, consequently, complex biological programs. In cancer, aberrant miRNA expression contributes to key oncogenic processes, including epithelial-mesenchymal transition, angiogenesis, immune evasion, and metastatic progression. Oncogenic miRNAs that lose tissue specificity during malignant transformation may represent promising therapeutic targets due to their limited expression in normal tissues. We performed a pan-cancer analysis combining tissue specificity profiling of miRNA expression in healthy tissues from the Genotype-Tissue Expression (GTEx) Project and tumor tissues from The Cancer Genome Atlas (TCGA) with differential expression analysis in tumor versus matched normal samples. Cross-cohort integration was used to identify miRNAs that both lose tissue specificity and display recurrent dysregulation. 

Malignant transformation was associated with a widespread loss of tissue-specific miRNA expression. Among these miRNAs, we identified a high-confidence cluster of nine oncomiRs: miR-105-5p, miR-1269a, miR-196a-5p, miR-9-5p, miR-96-5p, miR-210-3p, miR-301b-3p, miR-592, and miR-135b-5p which were significantly upregulated more frequently than downregulated across tumors. Functional enrichment analysis of experimentally validated targets indicated convergence on shared oncogenic pathways related to epithelial-mesenchymal transition, angiogenesis, hypoxia response, PI3K/AKT signaling, and immune modulation.


---

## Repository Structure

| File | Description |
|------|-------------|
| `load_data.R` | Data loading and processing workflow for the discovery and validation datasets |
| `functions.R` | auxillary functions |
| `cellchat.R` | **CellChat** pipeline for ligand–receptor interaction inference |
| `GAM.R` | Processing workflow for generative additive models implementation |
| `ESCAPE.R` | Workflow for functional enrichment analysis |
| `Diffusion_maps.R` | Workflow for diffusion maps and time trajectory analysis |
| `CNA.R` | InferCNA - copy number alteration inferring framework |
| `for_deconvolution.Rmd` | Script for data preparation for bulk-deconvolution analysis |
| `gse190504_CAF_rank_score.py` | Python script for CAF signature rank expression evaluation across bulk gse190504 data conditions |
| `gse190504_cellanneal_deconvolution.py` | Python script for cellanneal bulk-deconvolution analysis within bulk gse190504 data |

---


---

## Contact
Work email: amismailov@hse.ru
Personal email: neuro.promotion@gmail.com 

---
