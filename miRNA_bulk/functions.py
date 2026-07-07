import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import matplotlib.pyplot as plt
from scipy.stats import binomtest, mannwhitneyu, false_discovery_control
from zarr.core.sync import P
import auxillary_plots
import pathlib
import specificity_metrics
from constants import *


from sklearn.metrics import roc_auc_score, balanced_accuracy_score, roc_curve
from scipy.stats import mannwhitneyu
from cliffs_delta import cliffs_delta
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score



def argregate_gtex_meta(df_meta_in: pd.DataFrame) -> pd.DataFrame:
    df_meta = df_meta_in.copy(deep=True)
    df_meta['SMTSD'] = df_meta['SMTSD'].map(MAP_GTEX_META).fillna(df_meta['SMTSD'])
    df_meta = df_meta[df_meta['SMTSD'] != 'to_remove']
    return df_meta

def impute_tcga_meta(df_meta_in: pd.DataFrame) -> pd.DataFrame:
    # not reported imputation and filtering
    df_meta = df_meta_in.copy(deep=True)
    df_meta = df_meta[df_meta['Cancer_group'] != 'Not Reported']
    # 1. Словарь для обычных опухолей (где по названию рака понятен орган)
    map_tissue_organ = {
        'Glioma' : 'Brain',
        'Pheochromocytoma' : 'Adrenal Gland',
        'Adenocarcinoma' : 'Colon',
        'Myeloid leukemia' : 'Bone Marrow',
        'Germ cell tumor' : 'Testis',
        'B cell lymphoma' : 'Lymph nodes'
    }
    
    # 2. Словарь исключительно для здоровых тканей (Normal tissue), 
    # так как они распределены по разным проектам/органам
    map_normal_project_id = {
        'TARGET-AML': 'Bone Marrow',
        'TCGA-UCEC': 'Uterus',
        'CGCI-BLGSP': 'Lymph nodes',
        'CPTAC-3': None,  # Если этот проект нужно отфильтровать или обработать отдельно
        'TCGA-GBM': 'Brain'
    }
    
    # Базовая маска: ищем строки, где орган не указан
    mask_not_reported = df_meta['Tissue_Organ'] == 'Not Reported'
    
    # Специфичные маски для разделения логики
    mask_normal = mask_not_reported & (df_meta['Cancer_group'] == 'Normal tissue')
    mask_tumor = mask_not_reported & (df_meta['Cancer_group'] != 'Normal tissue')
    
    # ШАГ 1: Для обычных опухолей восстанавливаем орган по Cancer_group
    df_meta.loc[mask_tumor, 'Tissue_Organ'] = df_meta.loc[mask_tumor, 'Cancer_group'].map(map_tissue_organ)
    
    # ШАГ 2: Для Normal tissue смотрим на Project ID и берем орган оттуда
    df_meta.loc[mask_normal, 'Tissue_Organ'] = df_meta.loc[mask_normal, 'Project ID'].map(map_normal_project_id)
    
    # ШАГ 3: Опционально — удаляем строки, для которых маппинг вернул None (например, CPTAC-3)
    # или если встретился неизвестный проект/рак, которого не было в словарях
    df_meta = df_meta.dropna(subset=['Tissue_Organ'])
    # На случай если где-то осталась строка "Not Reported", которую не смогли замапить:
    df_meta = df_meta[df_meta['Tissue_Organ'] != 'Not Reported']

    map_tissue_organ_rename = {
        'Gallbladder' : 'Liver–Gallbladder',
        'Liver' : 'Liver–Gallbladder',
        'Soft tissues' : 'Soft Tissues',
        'Soft tissue' : 'Soft Tissues',
        'Bone marrow' : 'Bone Marrow',
        'Prostate gland' : 'Prostate',
        'Thyroid gland' : 'Thyroid',
        'Peritoneum' : 'Peritoneum',
        'Oral mucosa' : 'Mucosa',
        'Eye' : 'to_remove',
        'Thorax' : 'Thorax',
        'Thymus' : 'Thymus',
        'Adrenal gland' : 'Adrenal Gland'
    }

    df_meta['Tissue_Organ'] = df_meta['Tissue_Organ'].map(map_tissue_organ_rename).fillna(df_meta['Tissue_Organ'])
    df_meta = df_meta[df_meta['Tissue_Organ'] != 'to_remove']
    return df_meta

def standardize_name(name) -> str:
    if isinstance(name, str):
        return name.replace("_", "-").lower()
    return "".join(n.replace("_", "-").lower() for n in name)

def shared_mirs(path_tcga_counts=PATH_TCGA_COUNTS, path_gtex_counts=PATH_GTEX_COUNTS):
    tcga_counts = pd.read_parquet(path_tcga_counts)
    gtex_counts = pd.read_parquet(path_gtex_counts)
    tcga_counts['mature_name'] = tcga_counts['mature_name'].apply(standardize_name)
    gtex_counts['mature_name'] = gtex_counts['mature_name'].apply(standardize_name)
    shared_mirs = set(tcga_counts['mature_name']) & set(gtex_counts['mature_name'])
    return sorted(list(shared_mirs))

def _remove_metadata(df_in: pd.DataFrame) -> pd.DataFrame:
    df = df_in.copy(deep=True)
    meta_cols = [c for c in (META_TCGA_COLS + META_GTEX_COLS) if c in df.columns]
        
    df.drop(columns=meta_cols, inplace=True)
    return df

def _log2(df_in: pd.DataFrame) -> pd.DataFrame:
    df = df_in.copy(deep=True)

    # Сохраняем агрегатные колонки
    agg_cols = df[[c for c in [GTEX_AGGREGATE_GROUP, TCGA_AGGREGATE_GROUP] if c in df.columns]]

    # Удаляем их перед логарифмированием
    df = _remove_metadata(df)


    arr = np.log2(df.to_numpy(dtype=np.float64) + 1)
    df_log = pd.DataFrame(arr, index=df.index, columns=df.columns)

    # Возвращаем агрегатные колонки
    return pd.concat([df_log, agg_cols], axis=1)

def _tissue_expression_profile(
    df: pd.DataFrame,
    group_col: str,
    mir_cols: list[str],
    *,
    on_log: bool | None = None,
    log_mode: str | None = None,
) -> tuple[pd.DataFrame, pd.Series]:
    """
    Tissue/organ mean profile for specificity metrics.
    max_expr is always on linear mean-CPM scale (for min_CPM filtering).
    """
    on_log = SPECIFICITY_ON_LOG if on_log is None else on_log
    log_mode = SPECIFICITY_LOG_MODE if log_mode is None else log_mode

    profile_lin = df.groupby(group_col, dropna=False)[mir_cols].mean()
    max_expr = profile_lin.max(axis=0)

    if not on_log:
        return profile_lin, max_expr

    if log_mode == 'log_of_mean':
        profile = np.log2(profile_lin + 1.0)
    elif log_mode == 'mean_of_log':
        log_arr = np.log2(df[mir_cols].to_numpy(dtype=np.float64) + 1.0)
        log_df = pd.DataFrame(log_arr, index=df.index, columns=mir_cols)
        log_df[group_col] = df[group_col].values
        profile = log_df.groupby(group_col, dropna=False)[mir_cols].mean()
    else:
        raise ValueError("SPECIFICITY_LOG_MODE must be 'mean_of_log' or 'log_of_mean'")

    return profile, max_expr

def specificity_metric(x: np.ndarray) -> float:
    # TAU SPECIFICITY SCORE
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    x = np.clip(x, 0.0, None)

    n = x.size
    if n < 2:
        return np.nan

    mx = np.max(x)
    if mx <= 0 or not np.isfinite(mx):
        return 0.0   

    return float(np.sum(1.0 - x / mx) / (n - 1))

def get_top_tissues_expression(row, expression_df, threshold=0.85, tissue_col=None):
    """
    Function to find the top 3 tissues/cancer groups with the highest expression.
    """
    specificity_col = [c for c in row.index if SPECIFICITY_METRIC in c][0]
    
    if pd.isna(row[specificity_col]) or row[specificity_col] <= threshold:
        return None

    mirna = row["mature_name"]

    # 2. Адаптируемся под структуру датасета (TCGA или GTEx)
    if tissue_col and tissue_col in expression_df.columns:
        if mirna not in expression_df.columns:
            return None
        # Для TCGA считаем среднее по группам
        expr = expression_df.groupby(tissue_col)[mirna].mean()
    else:
        if mirna not in expression_df.columns:
            return None
        expr = expression_df[mirna]

    # 3. Считаем общую сумму
    total_expr = expr.sum()
    if total_expr == 0:
        return "Zero expression"

    # 4. Вычисляем доли для ВСЕХ классов
    all_fractions = expr / total_expr
    
    # --- УМНАЯ ФИЛЬТРАЦИЯ ---
    # Оставляем только те классы, где доля строго больше нуля,
    # сортируем по убыванию и берем МАКСИМУМ 3 штуки
    active_fractions = all_fractions[all_fractions > 0].sort_values(ascending=False).head(3)
    
    if active_fractions.empty:
        return "Zero expression"

    # Формируем красивую строку
    parts = [f"{cls} [{frac:.2f}]" for cls, frac in active_fractions.items()]
    return ", ".join(parts)


class Load_GTEx_Data:
    def __init__(self, path_gtex_counts, path_gtex_meta):
        self.path_gtex_counts = path_gtex_counts
        self.path_gtex_meta = path_gtex_meta
    
    def load_meta(self):
        meta = pd.read_parquet(self.path_gtex_meta)
        meta = meta.rename(columns={"SAMPID": "sample_id"})
        meta["sample_id"] = meta["sample_id"].astype(str)
        meta[GTEX_AGGREGATE_GROUP] = meta[GTEX_AGGREGATE_GROUP].astype(str)

        if GTEX_AGGREGATE_GROUP == 'SMTSD':
            meta = argregate_gtex_meta(meta)

        smts_n = meta.groupby(GTEX_AGGREGATE_GROUP).size()
        eligible_smts = smts_n[smts_n >= MIN_GTEX_SMTS_N].index
        meta_elig = meta[meta[GTEX_AGGREGATE_GROUP].isin(eligible_smts)].copy()
        print("Samples in meta:", len(meta), f"| after {GTEX_AGGREGATE_GROUP} ≥", MIN_GTEX_SMTS_N, ":", len(meta_elig))
        print(f"Number of eligible {GTEX_AGGREGATE_GROUP}:", len(eligible_smts))
        
        return meta_elig[['sample_id', GTEX_AGGREGATE_GROUP]]
    
    def load_counts(self, meta_elig=None):
        if meta_elig is None:
            meta_elig = self.load_meta()
        g = pd.read_parquet(self.path_gtex_counts)
        if "mature_name" not in g.columns:
            raise ValueError("mature_name column not found in GTEx counts")

        # standardize mature_name and retain only shared miRNAs
        g['mature_name'] = g['mature_name'].apply(standardize_name)
        g = g[g['mature_name'].isin(shared_mirs())]

        wide = g.set_index("mature_name").T

        wide.index = wide.index.astype(str)

        wide = wide.loc[wide.index.isin(meta_elig["sample_id"])]
        print("counts dataframe size (samples × miRNA):", wide.shape)

        # remove miRNAs with zero expression
        non_zero_mirnas = wide.sum(axis=0) > 0
        wide = wide.loc[:, non_zero_mirnas]

        smt_map = meta_elig.drop_duplicates("sample_id").set_index("sample_id")[GTEX_AGGREGATE_GROUP]
        smts = smt_map.reindex(wide.index)
        if smts.isna().any():
            wide = wide.loc[smts.notna()]
            smts = smts.loc[wide.index]
        print("with SMTS:", len(wide))
        return wide, smts
    
    def agregate(self, df: pd.DataFrame) -> pd.DataFrame:
        return df.groupby(GTEX_AGGREGATE_GROUP, dropna=False).mean()

    def CPM(self, wide=None, smts=None, log=False):
        if wide is None:
            wide, smts = self.load_counts()
        if smts is None:
            smts = self.load_meta()[GTEX_AGGREGATE_GROUP]
        
        counts = wide.astype(float)
        lib_size = counts.sum(axis=1)
        cpm = counts.div(lib_size, axis=0) * 1_000_000

        cpm[GTEX_AGGREGATE_GROUP] = smts.values

        if log:
            cpm = self.agregate(cpm)
            return _log2(cpm)
        
        return cpm

    def specificity_score(self, cpm_df, on_log: bool | None = None):
        on_log = SPECIFICITY_ON_LOG if on_log is None else on_log
        log_note = f", log2(CPM+1) [{SPECIFICITY_LOG_MODE}]" if on_log else ", linear CPM"
        print(f'Running {SPECIFICITY_METRIC} score, per {GTEX_AGGREGATE_GROUP}{log_note}...')
        mir_cols = [c for c in cpm_df.columns if c != GTEX_AGGREGATE_GROUP]
        profile, max_expr = _tissue_expression_profile(
            cpm_df, GTEX_AGGREGATE_GROUP, mir_cols, on_log=on_log
        )
        specificity_series = profile.apply(lambda col: specificity_metric(col.values), axis=0)
        specificity_series.loc[max_expr < min_CPM] = 0
        return (
            pd.DataFrame({"mature_name": specificity_series.index.astype(str), f"{SPECIFICITY_METRIC}_gtex": specificity_series.values})
            .sort_values(f"{SPECIFICITY_METRIC}_gtex", ascending=False, na_position="last")
            .reset_index(drop=True)
        )

    def run_pipeline(self):
        meta_elig = self.load_meta()
        wide, smts = self.load_counts(meta_elig)
        cpm = self.CPM(wide, smts)
        specificity_gtex_df = self.specificity_score(cpm)
        return cpm, specificity_gtex_df

    def plot_pie(self, title_prefix: str = "Samples", save_path: str = None):
        data = self.CPM(log=False)
        auxillary_plots.plot_pie_single(data, GTEX_AGGREGATE_GROUP, title_prefix, save_path)


class Load_TCGA_Data:
    def __init__(self, path_tcga_counts, path_tcga_meta):
        self.path_tcga_counts = path_tcga_counts
        self.path_tcga_meta = path_tcga_meta
        self.META_COLS = ["Tissue_Organ", "Sample Type", "Project ID", "Cancer_group"]

        self.MIN_TUMORS = 10
        self.MIN_NORMALS = 5

        self.NORMAL_SAMPLE_TYPES = {
            "Solid Tissue Normal",
            "Bone Marrow Normal",
            "Lymphoid Normal",
        }
        self.TUMOR_SAMPLE_TYPES = {
            "Primary Tumor",
            "Primary Blood Derived Cancer - Bone Marrow",
            "Primary Blood Derived Cancer - Peripheral Blood",
            'Recurrent Blood Derived Cancer - Bone Marrow',
            'Recurrent Blood Derived Cancer - Peripheral Blood',
            'Recurrent Tumor',
        }

        self.SAMPLE_TYPE = 'Sample Type'
        self.REMOVE_FROM_SAMPLE_TYPE = {
            'Blood Derived Cancer - Peripheral Blood, Post-treatment',
            'Blood Derived Cancer - Bone Marrow, Post-treatment'
        }

        # column with normal sample types: will be removed before comparison with GTEx
        self.TO_REMOVE_NORMAL_COLUMN = 'Sample Type'

        self.TISSUE_TYPE = 'Tissue_Organ'
        self.TISSUE_TYPE_TO_REMOVE = None #['Not Reported']

    def load_meta(self):
        meta = pd.read_parquet(self.path_tcga_meta)
        meta = impute_tcga_meta(meta)
        return meta
    
    def load_counts(self):
        counts = pd.read_parquet(self.path_tcga_counts)
        meta = self.load_meta()

        if "mature_name" not in counts.columns:
            raise ValueError("mature_name column not found in TCGA counts")

        # standardize mature_name and retain only shared miRNAs
        counts['mature_name'] = counts['mature_name'].apply(standardize_name)
        counts = counts[counts['mature_name'].isin(shared_mirs())]

        wide = counts.set_index("mature_name").T

        wide.index = wide.index.astype(str)
        wide.index.name = "sample_barcode"

        # метаданные строго по ключу sample_barcode (не по позиции строки!)
        m = meta.set_index(meta["sample_barcode"].astype(str))
        aligned = m.reindex(wide.index)

        for col in self.META_COLS:
            if col in aligned.columns:
                wide[col] = aligned[col].values

        _ann = [c for c in self.META_COLS if c in wide.columns]
        _mir = [c for c in wide.columns if c not in _ann]
        wide = wide[_ann + _mir]
        wide = wide[~wide[self.SAMPLE_TYPE].isin(self.REMOVE_FROM_SAMPLE_TYPE)]
        if self.TISSUE_TYPE_TO_REMOVE is not None:
            wide = wide[~wide[self.TISSUE_TYPE].isin(self.TISSUE_TYPE_TO_REMOVE)]
        return wide

    def agregate(self, df_in: pd.DataFrame, 
                aggregate_group: str = TCGA_AGGREGATE_GROUP, 
                remove_normals: bool = False) -> pd.DataFrame:

        df = df_in.copy(deep=True)
        meta_cols = [c for c in self.META_COLS if c in df.columns and c != aggregate_group]
        
        if remove_normals:
            df = df[~df[self.SAMPLE_TYPE].isin(self.NORMAL_SAMPLE_TYPES)]

        df.drop(columns=meta_cols, inplace=True)
        return df.groupby(aggregate_group, dropna=False).mean()

    def CPM(self, df_in=None, log=False):
        
        if df_in is None:
            wide = self.load_counts()
        else:
            wide = df_in.copy(deep=True)

        mir_cols = [c for c in wide.columns if c not in self.META_COLS]
        counts = wide[mir_cols].astype(float)
        lib_size = counts.sum(axis=1)
        cpm = counts.div(lib_size, axis=0) * 1_000_000
            
        cpm[self.META_COLS] = wide[self.META_COLS]

        vc = cpm[TCGA_AGGREGATE_GROUP].value_counts()
        keep = vc[vc > MIN_TCGA_Cancer_group_N].index
        cpm = cpm[cpm[TCGA_AGGREGATE_GROUP].isin(keep)]

        if log:
            cpm = self.agregate(cpm)
            return _log2(cpm)

        return cpm
    
    def specificity_score(self, lin_df_in, on_log: bool | None = None):
        lin_df = lin_df_in.copy(deep=True)
        on_log = SPECIFICITY_ON_LOG if on_log is None else on_log
        log_note = f", log2(CPM+1) [{SPECIFICITY_LOG_MODE}]" if on_log else ", linear CPM"

        print(f'Running {SPECIFICITY_METRIC} score, per {TCGA_AGGREGATE_GROUP}{log_note}...')
        meta_cols = [c for c in self.META_COLS if c in lin_df.columns]
        mir_cols = [c for c in lin_df.columns if c not in meta_cols and c != "sample_barcode"]
        profile, max_expr = _tissue_expression_profile(
            lin_df, TCGA_AGGREGATE_GROUP, mir_cols, on_log=on_log
        )
        specificity_series = profile.apply(lambda col: specificity_metric(col.values), axis=0)
        specificity_series.loc[max_expr < min_CPM] = 0
        specificity_tcga_df = (
            pd.DataFrame({"mature_name": specificity_series.index.astype(str), f"{SPECIFICITY_METRIC}_tcga": specificity_series.values})
            .sort_values(f"{SPECIFICITY_METRIC}_tcga", ascending=False, na_position="last")
            .reset_index(drop=True)
        )
        return specificity_tcga_df
    
    def run_pipeline(self):
        wide = self.load_counts()
        lin_df = self.CPM(wide)
        specificity_tcga_df = self.specificity_score(lin_df)
        return lin_df, specificity_tcga_df
    
    def plot_pie(self, col: str='Tissue_Organ', eligible_cols: list[str] = None, title_prefix: str = "TCGA samples", save_path: str = None):
        data = self.CPM(log=False)
        if eligible_cols is not None:
            data = data[data[col].isin(eligible_cols)]
        auxillary_plots.plot_pie_single(data, col, title_prefix, save_path)
    
    def tumor_vs_normal(self, df=None):
        # 1. Counts + annotation (wide: samples × miRNAs)
        if df is None:
            wide = self.load_counts()
        else:
            wide = df.copy(deep=True) # не мутируем!

        # index = sample_barcode (не полагаться на неявный порядок строк)
        if wide.index.name != "sample_barcode":
            if "sample_barcode" in wide.columns: # лучше проверить
                wide = wide.set_index("sample_barcode")
            wide.index.name = "sample_barcode"
        wide.index = wide.index.astype(str)

        present_meta = [c for c in self.META_COLS if c in wide.columns]
        if not present_meta:
            raise ValueError("No meta columns found in the data") # нет метаданных - гг

        mir_cols = [c for c in wide.columns if c not in present_meta]
        if not mir_cols:
            raise ValueError("No miRNA columns found in the data") # нет miRNA тем более - гг

        # 2. Metadata для DE: tumor / normal по Sample Type
        meta_de = self.load_meta().copy()

        def map_condition(sample_type) -> str | None:
            if pd.isna(sample_type):
                return None
            st = str(sample_type)
            if st in self.NORMAL_SAMPLE_TYPES:
                return "normal"
            if st in self.TUMOR_SAMPLE_TYPES:
                return "tumor"
            return None

        meta_de["sample_barcode"] = meta_de["sample_barcode"].astype(str)
        meta_de["condition"] = meta_de["Sample Type"].map(map_condition)
        meta_de = meta_de.dropna(subset=["condition"])

        # 3. Inner join counts ∩ meta: только образцы, которые есть в обоих источниках 
        count_barcodes = set(wide.index)
        meta_barcodes = set(meta_de["sample_barcode"])
        shared = sorted(count_barcodes & meta_barcodes)
        missing_in_counts = meta_barcodes - count_barcodes
        if missing_in_counts:
            raise ValueError(
                f"{len(missing_in_counts)} tumor/normal samples in meta "
                f"are missing from counts (inner join failed)."
            )
        if not shared:
            raise ValueError("No shared samples between counts and tumor/normal metadata.")

        meta_de = meta_de[meta_de["sample_barcode"].isin(shared)]
        wide = wide.loc[shared]

        # 4. Integer count matrix для DESeq2 
        arr = wide[mir_cols].to_numpy(dtype=float) # тут важно что бы wide был матрицей каунтов, иначе все сломается
        arr = np.maximum(arr, 0) # убираем отрицательные значения на всякий, но их и не должно быть
        if not np.allclose(arr, np.round(arr), rtol=0, atol=1e-6): 
            arr = np.round(arr)
        counts_all = pd.DataFrame(
            arr.astype(np.int64), # приводим к виду int явно, чтобы DESeq2 не ругался
            index=pd.Index(shared, name="sample_barcode"),
            columns=mir_cols,
        )
        # без этой колонки никак
        if "Project ID" not in meta_de.columns:
            raise ValueError("Project ID column not found in the data")

        # 5. Берем проекты с достаточным числом tumor и normal
        proj_tab = (
            meta_de.groupby(["Project ID", "condition"], observed=False)
            .size()
            .unstack(fill_value=0)
        )
        for c in ("normal", "tumor"):
            if c not in proj_tab.columns:
                proj_tab[c] = 0

        eligible_projects = proj_tab[
            (proj_tab["tumor"] >= self.MIN_TUMORS) & (proj_tab["normal"] >= self.MIN_NORMALS)
        ].index.tolist()

        if not eligible_projects:
            raise RuntimeError(
                f"No projects with >={self.MIN_TUMORS} tumors and >={self.MIN_NORMALS} normals."
            )

        # 6. DESeq2 tumor vs normal отдельно по каждому eligible-проекту (ну там где есть достаточное число tumor и normal)
        de_parts: list[pd.DataFrame] = []

        for project_id in eligible_projects:
            m = meta_de[meta_de["Project ID"] == project_id]
            ids_normal = m.loc[m["condition"] == "normal", "sample_barcode"].drop_duplicates().tolist()
            ids_tumor = m.loc[m["condition"] == "tumor", "sample_barcode"].drop_duplicates().tolist()

            if not ids_normal or not ids_tumor:
                continue

            sample_order = ids_normal + ids_tumor
            counts_sub = counts_all.loc[sample_order]

            meta_f = pd.DataFrame(
                {"condition": ["normal"] * len(ids_normal) + ["tumor"] * len(ids_tumor)},
                index=sample_order,
            )
            meta_f["condition"] = pd.Categorical(
                meta_f["condition"], categories=["normal", "tumor"]
            )

            # miRNA с counts >= 5 хотя бы в 3 образцах проекта
            mat = counts_sub.to_numpy().T
            keep = (mat >= 5).sum(axis=1) >= 3
            counts_f = counts_sub.loc[:, keep]

            if counts_f.shape[1] == 0:
                continue

            try:
                dds = DeseqDataSet(
                    counts=counts_f,
                    metadata=meta_f,
                    design="~condition",
                    refit_cooks=True,
                    quiet=True,
                )
                dds.deseq2()
                stat = DeseqStats(
                    dds,
                    contrast=["condition", "tumor", "normal"],
                    inference=DefaultInference(n_cpus=1),
                    quiet=True,
                )
                stat.summary()
            except Exception:
                continue

            res = stat.results_df.reset_index()
            if res.columns[0] != "mature_name":
                res = res.rename(columns={res.columns[0]: "mature_name"})

            res["project_id"] = project_id
            res["n_normal"] = len(ids_normal)
            res["n_tumor"] = len(ids_tumor)
            de_parts.append(res)

        if not de_parts:
            raise RuntimeError( # не дай бог
                "No successful DE results — check meta, counts, and eligible projects."
            )

        de_all = pd.concat(de_parts, ignore_index=True)
        return de_all, proj_tab, eligible_projects

    


def binom_pval(n_up: int, n_down: int) -> tuple[float, float, float]:

    n = int(n_up) + int(n_down)
    if n == 0:
        return (float("nan"), float("nan"), float("nan"))
    k = int(n_up)
    p2 = float(binomtest(k, n, p=0.5, alternative="two-sided").pvalue)
    p_up = float(binomtest(k, n, p=0.5, alternative="greater").pvalue)
    p_dn = float(binomtest(k, n, p=0.5, alternative="less").pvalue)
    return (p2, p_up, p_dn)




def compare_gtex_tcga_profiles( 
    alternative: str = "two-sided",
    fdr_method: str = "bh",
) -> pd.DataFrame:
    """
    Compare pan-tissue mean CPM profiles: GTEx (normal tissues) vs TCGA (tumor organs).

    Per miRNA: tissue/organ mean linear CPM → log2(CPM+1) → Mann–Whitney U on log values.
    log2FC = median(log2 TCGA) − median(log2 GTEx). direction requires |log2FC| > LFC_THR
    and padj_tcga_gtex < PADJ_THR.
    """
    if alternative not in {"two-sided", "less", "greater"}:
        raise ValueError("alternative must be 'two-sided', 'less', or 'greater'")

    gtex_loader = Load_GTEx_Data(PATH_GTEX_COUNTS, PATH_GTEX_META)
    tcga_loader = Load_TCGA_Data(PATH_TCGA_COUNTS, PATH_TCGA_META)

    gtex_cpm = gtex_loader.CPM(log=False)
    tcga_cpm = tcga_loader.CPM(log=False)

    gtex = gtex_loader.agregate(gtex_cpm)
    tcga = tcga_loader.agregate(tcga_cpm, remove_normals=True)

    if not gtex.columns.equals(tcga.columns):
        shared = [c for c in gtex.columns if c in tcga.columns]
        if not shared:
            raise ValueError("No shared miRNA columns between GTEx and TCGA profiles")
        gtex = gtex[shared]
        tcga = tcga[shared]

    mir_cols = list(gtex.columns)
    rows: list[dict] = []

    for mir in mir_cols:
        x = gtex[mir].to_numpy(dtype=float)
        y = tcga[mir].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]

        med_gtex = float(np.median(x)) if x.size else np.nan
        med_tcga = float(np.median(y)) if y.size else np.nan

        x_log = np.log2(x + 1.0)
        y_log = np.log2(y + 1.0)
        med_log_gtex = float(np.median(x_log)) if x_log.size else np.nan
        med_log_tcga = float(np.median(y_log)) if y_log.size else np.nan
        log2_fc = med_log_tcga - med_log_gtex if np.isfinite(med_log_gtex) and np.isfinite(med_log_tcga) else np.nan

        if x_log.size == 0 or y_log.size == 0:
            stat, pval = np.nan, np.nan
        else:
            res = mannwhitneyu(x_log, y_log, alternative=alternative)
            stat, pval = float(res.statistic), float(res.pvalue)

        rows.append(
            {
                "mature_name": mir,
                "median_cpm_gtex": med_gtex,
                "median_cpm_tcga": med_tcga,
                "mean_cpm_gtex": float(np.mean(x)) if x.size else np.nan,
                "mean_cpm_tcga": float(np.mean(y)) if y.size else np.nan,
                "log2FC_median_tcga_gtex": log2_fc,
                "mw_statistic": stat,
                "pvalue": pval,
            }
        )

    out = pd.DataFrame(rows)
    valid = out["pvalue"].notna()
    out["padj_tcga_gtex"] = np.nan
    if valid.any():
        out.loc[valid, "padj_tcga_gtex"] = false_discovery_control(
            out.loc[valid, "pvalue"].to_numpy(dtype=float),
            method=fdr_method,
        )

    def _direction(row) -> str:
        fc = row["log2FC_median_tcga_gtex"]
        padj = row["padj_tcga_gtex"]
        if pd.isna(fc) or pd.isna(padj):
            return "NA"
        if fc > LFC_THR_CPM and padj < PADJ_THR:
            return "TCGA_UP"
        if fc < -LFC_THR_CPM and padj < PADJ_THR:
            return "GTEx_UP"
        return "NS"

    out["direction"] = out.apply(_direction, axis=1)
    return out.sort_values("pvalue", na_position="last").reset_index(drop=True)


def _strip_ensembl_version(gene_id) -> str:
    return str(gene_id).split(".")[0]


def _gtex_sample_key(sample_id: str) -> str:
    parts = str(sample_id).split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else str(sample_id)


def _map_gtex_mirna_to_mrna_samples(mir_samples, rna_samples) -> dict[str, str]:
    """GTEx miRNA and mRNA use different aliquot suffixes; match on GTEX-XXX-YYYY."""
    rna_by_key: dict[str, str] = {}
    for s in rna_samples:
        k = _gtex_sample_key(s)
        if k not in rna_by_key:
            rna_by_key[k] = s
    return {s: rna_by_key[_gtex_sample_key(s)] for s in mir_samples if _gtex_sample_key(s) in rna_by_key}


def _read_gtex_mrna_sample_columns(path=PATH_GTEX_MRNA) -> list[str]:
    import gzip
    with gzip.open(path, "rt") as f:
        f.readline()
        f.readline()
        header = f.readline().strip().split("\t")
    return header[2:]


def _read_tcga_mrna_sample_columns(path=PATH_TCGA_MRNA) -> list[str]:
    return pd.read_csv(path, nrows=0).columns.tolist()[1:]


def _load_gtex_mrna(path=PATH_GTEX_MRNA, columns: list[str] | None = None) -> pd.DataFrame:
    usecols = None
    if columns is not None:
        usecols = ["Name", "Description", *columns]
    df = pd.read_csv(path, sep="\t", skiprows=2, low_memory=False, usecols=usecols)
    df["Name"] = df["Name"].map(_strip_ensembl_version)
    df = df.set_index("Name")
    df.index = df.index.astype(str)
    df = df.drop(columns=["Description"], errors="ignore")
    return df


def _load_tcga_mrna(path=PATH_TCGA_MRNA, columns: list[str] | None = None) -> pd.DataFrame:
    usecols = None
    if columns is not None:
        usecols = ["gene_id", *columns]
    df = pd.read_csv(path, low_memory=False, usecols=usecols)
    df["gene_id"] = df["gene_id"].map(_strip_ensembl_version)
    df = df.set_index("gene_id")
    df.index = df.index.astype(str)
    return df


def _counts_to_log2_cpm(
    df: pd.DataFrame,
    lib_size: pd.Series | None = None,
) -> pd.DataFrame:
    counts = df.astype(np.float64)
    if lib_size is None:
        lib = counts.sum(axis=1)
    else:
        lib = lib_size.reindex(counts.index)
    lib = lib.replace(0, np.nan)
    cpm = counts.div(lib, axis=0) * 1_000_000
    return np.log2(cpm + 1.0)


def _zscore_matrix(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Column-wise z-score; zero-variance columns marked in valid mask."""
    mu = np.nanmean(arr, axis=0)
    sd = np.nanstd(arr, axis=0, ddof=1)
    valid = sd > 0
    z = np.zeros_like(arr, dtype=np.float64)
    z[:, valid] = (arr[:, valid] - mu[valid]) / sd[valid]
    return z, valid


def _top_anticorr_genes(
    mir_mat: np.ndarray,
    rna_mat: np.ndarray,
    mir_names: list[str],
    gene_names: list[str],
    top_k: int,
) -> pd.DataFrame:
    """Top-K mRNA genes with the lowest Pearson r (strongest anti-correlation)."""
    mir_z, mir_ok = _zscore_matrix(mir_mat)
    rna_z, rna_ok = _zscore_matrix(rna_mat)
    n = mir_mat.shape[0] - 1
    if n <= 0:
        raise ValueError("Need at least 2 shared samples for correlation.")
    if top_k <= 0:
        raise ValueError("top_k must be a positive integer.")

    rna_ok_idx = np.where(rna_ok)[0]
    if rna_ok_idx.size == 0:
        return pd.DataFrame(columns=["mature_name", "gene_id", "pearson_r", "rank"])

    rna_z_ok = rna_z[:, rna_ok_idx]
    genes_ok = [gene_names[i] for i in rna_ok_idx]
    corr = (mir_z[:, mir_ok].T @ rna_z_ok) / n
    mir_ok_names = [mir_names[i] for i, ok in enumerate(mir_ok) if ok]
    k = min(top_k, len(genes_ok))

    rows: list[dict] = []
    for mir, row in zip(mir_ok_names, corr):
        finite = np.isfinite(row)
        if not finite.any():
            continue
        row_f = row.copy()
        row_f[~finite] = np.inf
        idx = np.argpartition(row_f, k - 1)[:k]
        idx = idx[np.argsort(row_f[idx])]
        for rank, gene_i in enumerate(idx, start=1):
            rows.append(
                {
                    "mature_name": mir,
                    "gene_id": genes_ok[gene_i],
                    "pearson_r": float(row[gene_i]),
                    "rank": rank,
                }
            )

    if not rows:
        return pd.DataFrame(columns=["mature_name", "gene_id", "pearson_r", "rank"])
    return pd.DataFrame(rows)


def anti_correlation_analysis(
    mirs: list[str],
    top_k: int = 100,
) -> pd.DataFrame:
    """
    For each miRNA, return top-K mRNA genes with the strongest anti-correlation
    across GTEx + TCGA samples (lowest Pearson r).

    Expression: log2(CPM + 1). Returns a DataFrame with columns:
    mature_name, gene_id, pearson_r, rank (1 = strongest anti-correlation).
    """
    if not mirs:
        return pd.DataFrame(columns=["mature_name", "gene_id", "pearson_r", "rank"])

    mirs = [standardize_name(m) for m in mirs]

    gtex_loader = Load_GTEx_Data(PATH_GTEX_COUNTS, PATH_GTEX_META)
    tcga_loader = Load_TCGA_Data(PATH_TCGA_COUNTS, PATH_TCGA_META)

    gtex_mir = gtex_loader.load_counts()[0]
    tcga_mir = _remove_metadata(tcga_loader.load_counts())

    for mir in mirs:
        in_gtex = mir in gtex_mir.columns
        in_tcga = mir in tcga_mir.columns
        if not in_gtex and not in_tcga:
            raise ValueError(f"miRNA {mir} not found in GTEx or TCGA")

    gtex_map = _map_gtex_mirna_to_mrna_samples(
        gtex_mir.index, _read_gtex_mrna_sample_columns()
    )
    gtex_mir_samples = list(gtex_map.keys())
    gtex_rna_samples = [gtex_map[s] for s in gtex_mir_samples]
    tcga_shared_samples = tcga_mir.index.intersection(_read_tcga_mrna_sample_columns())

    df_rna_gtex = _load_gtex_mrna(columns=gtex_rna_samples)
    df_rna_tcga = _load_tcga_mrna(columns=list(tcga_shared_samples))

    shared_genes = df_rna_gtex.index.intersection(df_rna_tcga.index)
    df_rna_gtex = df_rna_gtex.loc[shared_genes]
    df_rna_tcga = df_rna_tcga.loc[shared_genes]

    gtex_mir_sub = gtex_mir.loc[gtex_mir_samples, [m for m in mirs if m in gtex_mir.columns]]
    gtex_rna_sub = df_rna_gtex[gtex_rna_samples].T
    gtex_rna_sub.index = gtex_mir_samples

    tcga_mir_sub = tcga_mir.loc[tcga_shared_samples, [m for m in mirs if m in tcga_mir.columns]]
    tcga_rna_sub = df_rna_tcga[list(tcga_shared_samples)].T
    tcga_rna_sub.index = tcga_shared_samples

    df_mir = pd.concat([gtex_mir_sub, tcga_mir_sub], axis=0)
    df_rna = pd.concat([gtex_rna_sub, tcga_rna_sub], axis=0)

    missing_mirs = [m for m in mirs if m not in df_mir.columns]
    if missing_mirs:
        raise ValueError(f"miRNAs missing after sample alignment: {missing_mirs}")

    df_mir = df_mir[mirs]
    if not df_mir.index.equals(df_rna.index):
        shared_samples = df_mir.index.intersection(df_rna.index)
        df_mir = df_mir.loc[shared_samples]
        df_rna = df_rna.loc[shared_samples]

    if len(df_mir) < 2:
        raise ValueError("Fewer than 2 aligned samples for correlation.")

    # CPM library size from full count matrices, not from selected miRNA columns only
    gtex_lib = gtex_mir.loc[gtex_mir_samples].sum(axis=1)
    tcga_lib = tcga_mir.loc[tcga_shared_samples].sum(axis=1)
    mir_lib = pd.concat([gtex_lib, tcga_lib]).reindex(df_mir.index)

    df_mir = _counts_to_log2_cpm(df_mir, lib_size=mir_lib)
    df_rna = _counts_to_log2_cpm(df_rna)
    df_mir = df_mir.fillna(0.0)
    df_rna = df_rna.fillna(0.0)

    return _top_anticorr_genes(
        df_mir.to_numpy(dtype=np.float64),
        df_rna.to_numpy(dtype=np.float64),
        mir_names=list(df_mir.columns),
        gene_names=list(df_rna.columns),
        top_k=top_k,
    )


def classification_statistics(df_in):
    df = df_in[df_in["Dysregulated_class"].isin(["DOWNREGULATED", "UPREGULATED"])].copy()
    y = (df["Dysregulated_class"] == "UPREGULATED").astype(int)
    scores = -df["Delta_Tau"]
    
    # Расчеты
    auc = roc_auc_score(y, scores)
    
    up = df.loc[y == 1, "Delta_Tau"]
    down = df.loc[y == 0, "Delta_Tau"]
    stat, p_val = mannwhitneyu(up, down, alternative='two-sided')
    delta, _ = cliffs_delta(up, down)
    
    # Оптимальный порог
    fpr, tpr, thresholds = roc_curve(y, scores)
    j_scores = tpr - fpr
    optimal_threshold = thresholds[np.argmax(j_scores)]
    
    # ПРЕДСКАЗАНИЯ по оптимальному порогу (вместо CV для малых выборок)
    # Это гораздо честнее для 117 объектов
    y_pred = (scores >= optimal_threshold).astype(int)
    bal_acc = balanced_accuracy_score(y, y_pred)
    
    print("--- Classification Statistics Report ---")
    print(f"ROC-AUC:               {auc:.3f}")
    print(f"Balanced Accuracy:     {bal_acc:.3f} (at optimal threshold)")
    print(f"Optimal Threshold:     {optimal_threshold:.3f}")
    print(f"Mann-Whitney U p-value: {p_val:.2e}")
    print(f"Cliff's Delta:         {delta:.3f}")
    print("----------------------------------------")
    return auc