import pathlib
PROJ_ROOT = pathlib.Path("/Users/neuropromotion/Desktop/miRNA/bulk_review")

PATH_GTEX_COUNTS = PROJ_ROOT / "data/processed/gtex_mature_counts_filtered.parquet"
PATH_GTEX_META = PROJ_ROOT / "data/processed/gtex_mature_metadata.parquet"

PATH_TCGA_META = PROJ_ROOT / "data/processed/tcga_mature_metadata.parquet"
PATH_TCGA_COUNTS = PROJ_ROOT / "data/processed/tcga_mature_counts_filtered.parquet"

PATH_TCGA_MRNA = PROJ_ROOT / 'data/TCGA_genes_counts.csv'
PATH_GTEX_MRNA = PROJ_ROOT / 'data/GTEX_genes_counts.gct.gz'

FIG_OUT = PROJ_ROOT / "pipeline_python" / "figures"
FIG_OUT.mkdir(parents=True, exist_ok=True)

TAB_OUT = PROJ_ROOT / "pipeline_python" / "tables"
TAB_OUT.mkdir(parents=True, exist_ok=True)

META_TCGA_COLS = ["Tissue_Organ", "Sample Type", "Project ID", "Cancer_group"]
META_GTEX_COLS = ["SMTSD", "SMTS"]

MIN_GTEX_SMTS_N = 10
MIN_TCGA_Cancer_group_N = 10

GTEX_AGGREGATE_GROUP = 'SMTSD' # 'SMTSD' | 'SMTS'
TCGA_AGGREGATE_GROUP = 'Tissue_Organ' #'Tissue_Organ' 'Cancer_group'

SPECIFICITY_METRIC = 'Tau'
SPECIFICITY_THR = 0.8  

SPECIFICITY_ON_LOG = True
SPECIFICITY_LOG_MODE = 'log_of_mean' 

BINOM_P_THR = 0.05 
PADJ_THR = 0.05 # threshold for significant DE
LFC_THR = 1.0 # threshold for significant DE
LFC_THR_CPM = 1

min_CPM = -1 # shutted down 
FIG_MODE = 'pdf'

TOP_N = 100 # shutted down 
OTHER = "Other types"
MISSING = "<missing>"


FINAL_CHOICES = [
        'concordant specific miR',
        'cancer-tissue-specific miR',
        'normal-tissue-specific miR',
        'concordant non-specific miR'
    ]

FINAL_COLORS = {
    'concordant specific miR': '#2ca02c',
    'cancer-tissue-specific miR': '#9467bd',
    'normal-tissue-specific miR': '#B8860B',
    'concordant non-specific miR': '#7f7f7f'
}


STN_CLASS_COLORS = {
    "NOT_SIGNIFICANT": "gray",
    "UPREGULATED": "firebrick",
    "DOWNREGULATED": "navy",
}

SPECIFICITY_AXIS_PAD = 0.05

MAP_GTEX_META = {
        'Cells - EBV-transformed lymphocytes' : 'to_remove',
        'Cells - Cultured fibroblasts' : 'to_remove',
        'Kidney - Cortex' : 'Kidney',
        'Kidney - Medulla' : 'Kidney',
        'Esophagus - Mucosa' : 'Mucosa',
        'Esophagus - Muscularis' : 'Esophagus',
        'Esophagus - Gastroesophageal Junction' : 'Esophagus',

        'Small Intestine - Terminal Ileum' : 'Colon–Intestine',
        'Colon - Sigmoid' : 'Colon–Intestine',
        'Colon - Transverse' : 'Colon–Intestine',
        
        'Brain - Cortex' : 'Brain',
        'Brain - Caudate (basal ganglia)' : 'Brain',
        'Brain - Putamen (basal ganglia)' : 'Brain',
        'Brain - Substantia nigra' : 'Brain',
        'Brain - Hippocampus' : 'Brain',
        'Brain - Spinal cord (cervical c-1)' : 'Brain',
        'Brain - Cerebellum' : 'Brain',
        'Brain - Nucleus accumbens (basal ganglia)' : 'Brain',
        'Brain - Hypothalamus' : 'Brain',
        'Brain - Frontal Cortex (BA9)' : 'Brain',
        'Brain - Anterior cingulate cortex (BA24)' : 'Brain',
        'Brain - Cerebellar Hemisphere' : 'Brain',
        'Brain - Amygdala' : 'Brain',
        'Pituitary' : 'Brain',

        'Heart - Left Ventricle' : 'Heart',
        'Heart - Atrial Appendage' : 'Heart',

        'Artery - Tibial' : 'Soft Tissues',
        'Artery - Coronary' : 'Soft Tissues',
        'Artery - Aorta' : 'Soft Tissues',
        'Adipose - Subcutaneous' : 'Soft Tissues',
        'Nerve - Tibial' : 'Soft Tissues',
        'Muscle - Skeletal' : 'Soft Tissues',


        'Cervix - Ectocervix' : 'Uterus–Cervix',
        'Cervix - Endocervix' : 'Uterus–Cervix',
        'Uterus' : 'Uterus–Cervix',

        'Skin - Not Sun Exposed (Suprapubic)' : 'Skin',
        'Skin - Sun Exposed (Lower leg)' : 'Skin',

        'Adipose - Visceral (Omentum)' : 'Omentum',
        'Breast - Mammary Tissue' : 'Breast',

        'Minor Salivary Gland' : 'Salivary Gland',

        'Vagina' : 'Vagina'
    }