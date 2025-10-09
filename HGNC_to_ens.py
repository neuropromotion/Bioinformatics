from mygene import MyGeneInfo
import pandas as pd
from tqdm import tqdm

def find_ens(list_of_genes):
    mg = MyGeneInfo()
    
    # 1) Запрос по списку символов
    query = mg.querymany(
        list_of_genes,
        scopes='symbol',
        fields='ensembl.gene',
        species='human',
        as_dataframe=False
    )
    results = {}
    for entry in tqdm(query):
        sym = entry.get('query')
        ensembl = entry.get('ensembl')
    
        gene_id = None
        # ensembl может быть:
        # 1) словарём {'gene': 'ENSG...', ...}
        # 2) списком словарей [{'gene':'ENSG...', ...}, {...}]
        if isinstance(ensembl, dict):
            gene_id = ensembl.get('gene')
        elif isinstance(ensembl, list) and len(ensembl) > 0:
            gene_id = ensembl[0].get('gene')
        
        # если не нашли, оставляем None
        results[sym] = gene_id
    return results 
