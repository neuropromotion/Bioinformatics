# Seminar: Glioblastoma Single-Cell RNA-Seq Cell Classifier

Materials for a practical seminar dedicated to single-cell transcriptomics (scRNA-seq) data analysis and building a classifier for glioblastoma tumor and microenvironment cells.

---

## Structure

* `seminar.ipynb` — Main Jupyter Notebook.
* `chromosome_means.py` — Utility script for aggregating gene expression into chromosomal means (used for copy number variation / CNV estimation).
* `genes_chr_mapping.csv` — Mapping table linking genes to their respective chromosomes.
* `singlecell.yml` — Conda configuration file containing all necessary environment dependencies.
* `5v1/` — Target directory for the first dataset.
* `3v3/` — Target directory for the second dataset.

---

## Environment Setup

To run the seminar code smoothly, it is recommended to create a dedicated Conda environment using the provided configuration file:

```bash
# Create the environment
conda env create -f singlecell.yml

# Activate the environment
conda activate <environment_name>
```
## Data download links:
5v1 scRNA-seq dataset: https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-5-v-1-targeted-neuroscience-panel-1-standard-4-0-0
3v3 scRNA-seq dataset: https://www.10xgenomics.com/datasets/human-glioblastoma-multiforme-3-v-3-whole-transcriptome-analysis-3-standard-4-0-0
