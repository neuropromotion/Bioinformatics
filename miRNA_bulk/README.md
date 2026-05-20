# Loss of tissue specificity and recurrent pan-cancer activation define a conserved oncogenic microRNA program

Overview

MicroRNAs are post-transcriptional regulators capable of coordinately modulating large gene networks and, consequently, complex biological programs. In cancer, aberrant miRNA expression contributes to key oncogenic processes, including epithelial-mesenchymal transition, angiogenesis, immune evasion, and metastatic progression. Oncogenic miRNAs that lose tissue specificity during malignant transformation may represent promising therapeutic targets due to their limited expression in normal tissues. We performed a pan-cancer analysis combining tissue specificity profiling of miRNA expression in healthy tissues from the Genotype-Tissue Expression (GTEx) Project and tumor tissues from The Cancer Genome Atlas (TCGA) with differential expression analysis in tumor versus matched normal samples. Cross-cohort integration was used to identify miRNAs that both lose tissue specificity and display recurrent dysregulation. 

Malignant transformation was associated with a widespread loss of tissue-specific miRNA expression. Among these miRNAs, we identified a high-confidence cluster of nine oncomiRs: miR-105-5p, miR-1269a, miR-196a-5p, miR-9-5p, miR-96-5p, miR-210-3p, miR-301b-3p, miR-592, and miR-135b-5p which were significantly upregulated more frequently than downregulated across tumors. Functional enrichment analysis of experimentally validated targets indicated convergence on shared oncogenic pathways related to epithelial-mesenchymal transition, angiogenesis, hypoxia response, PI3K/AKT signaling, and immune modulation.


---

## Repository Structure

| File | Description |
|------|-------------|
| `TCGA.ipynb` | Notebook with TCGA dataset processing: Tau specificity assessment, DE analysis (cancer-vs-normal) and finally binomial test implementation |
| `GTEx.ipynb` | Notebook with Tau specificity assessment within GTEx dataset (normal tissues) |
| `TCGA_GTEx.ipynb` | Integration specificity (Tau metric) and reccurent dysregulation |
| `Interactive_2D_results_map.html` | Interactive map with all integrated results | 

---


---

## Contact
Work email: amismailov@hse.ru
Personal email: neuro.promotion@gmail.com 

---
