import time
import requests
import pandas as pd
from io import StringIO

STRING_API = "https://version-12-0.string-db.org/api"  

STRING_CATEGORIES = {
    "Process": "GO: Biological Process",
    "Function": "GO: Molecular Function",
    "Component": "GO: Cellular Component",
    "KEGG": "KEGG",
    "RCTM": "Reactome",
    "WikiPathways": "WikiPathways",
}

def string_enrichment(
    genes,
    species=9606,
    categories=("Process", "KEGG", "RCTM"),   # GO BP + KEGG + Reactome, можно и другие добавить, но я подумал что этого достаточно..
    fdr_thr=0.05,
    caller_identity="mirna_pipeline",
):
    genes = list(dict.fromkeys(g for g in genes if pd.notna(g)))   
    if len(genes) < 3:
        raise ValueError("STRING enrichment требует хотя бы ~3 гена.")

    resp = requests.post(
        f"{STRING_API}/tsv/enrichment",
        data={
            "identifiers": "\r".join(genes),    
            "species": species,
            "caller_identity": caller_identity,
        },
        timeout=60,
    )
    resp.raise_for_status()

    df = pd.read_csv(StringIO(resp.text), sep="\t")
    if df.empty:
        return df

    if categories:
        df = df[df["category"].isin(categories)]
    if fdr_thr is not None:
        df = df[df["fdr"] <= fdr_thr]

    df["source"] = df["category"].map(STRING_CATEGORIES).fillna(df["category"])
    # кратность обогащения: доля в списке / доля в фоне
    df["fold_enrichment"] = (
        (df["number_of_genes"] / len(genes))
        / (df["number_of_genes_in_background"] / df["number_of_genes_in_background"].max())
    )
    return df.sort_values(["category", "fdr"]).reset_index(drop=True)


def enrich_targets_per_mirna(target_lists: dict, pause=1.0, **kwargs):
    """target_lists: {'hsa-mir-421': [gene1, gene2, ...], ...}"""
    out = []
    for mir, genes in target_lists.items():
        try:
            df = string_enrichment(genes, **kwargs)
        except Exception as e:
            print(f"{mir}: пропущено ({e})")
            continue
        if not df.empty:
            df.insert(0, "miRNA", mir)
            out.append(df)
        time.sleep(pause)   # чтоб не перегрузить апишку
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()