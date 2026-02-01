from mygene import MyGeneInfo
import pandas as pd
from tqdm import tqdm
import requests

 def find_ens(self, genes, sleep=0.2, verbose=False):
        """
        Robust gene symbol / alias → ENSG mapping.
        Uses mygene.info + Ensembl REST fallback.
        """
    
        mg = MyGeneInfo()
        results = {}
    
        def extract_ensg(hit):
            ens = hit.get("ensembl")
            if isinstance(ens, dict):
                g = ens.get("gene")
                if isinstance(g, str) and g.startswith("ENSG"):
                    return g
            elif isinstance(ens, list):
                for e in ens:
                    g = e.get("gene")
                    if isinstance(g, str) and g.startswith("ENSG"):
                        return g
            return None
    
        def ensembl_lookup(symbol):
            url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{symbol}"
            headers = {"Content-Type": "application/json"}
            r = requests.get(url, headers=headers, timeout=10)
            if not r.ok:
                return None
            for entry in r.json():
                if entry.get("type") == "gene" and entry.get("id", "").startswith("ENSG"):
                    return entry["id"]
            return None
    
        for gene in genes:
            ensg = None
    
            # --- 1. mygene.info ---
            try:
                res = mg.query(
                    gene,
                    scopes="_all",
                    species="human",
                    fields="ensembl",
                    size=20
                )
                for hit in res.get("hits", []):
                    ensg = extract_ensg(hit)
                    if ensg:
                        break
            except Exception:
                pass
    
            # --- 2. Ensembl fallback ---
            if ensg is None:
                ensg = ensembl_lookup(gene)
    
            if verbose:
                print(f"{gene} → {ensg}")
    
            results[gene] = ensg
            time.sleep(sleep)
    
        return results
