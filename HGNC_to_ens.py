from mygene import MyGeneInfo
import pandas as pd
from tqdm import tqdm

def find_ens(list_of_genes, use_synonyms=True, verbose=False):
    """
    Конвертирует символы генов в Ensembl ID.
    При неудаче пытается найти через синонимы.
    
    Parameters:
    -----------
    list_of_genes : list
        Список символов генов
    use_synonyms : bool
        Использовать ли поиск по синонимам при неудаче (по умолчанию True)
    verbose : bool
        Выводить ли детали поиска по синонимам
        
    Returns:
    --------
    dict : {gene_symbol: ensembl_id or None}
    """
    def _search_by_synonyms(gene_symbol, mg, verbose=False):
        """
        Вспомогательная функция: ищет Ensembl ID через синонимы гена
        
        Parameters:
        -----------
        gene_symbol : str
            Символ гена, который не был найден
        mg : MyGeneInfo
            Экземпляр MyGeneInfo
        verbose : bool
            Выводить ли детальную информацию
            
        Returns:
        --------
        str or None : Ensembl gene ID или None
        """
        # Сначала получаем синонимы
        try:
            synonym_query = mg.query(
                gene_symbol,
                species='human',
                fields='symbol,alias,other_names,ensembl.gene',
                size=1
            )
            
            if not synonym_query or 'hits' not in synonym_query:
                return None
            
            hits = synonym_query['hits']
            if not hits:
                return None
            
            hit = hits[0]
            
            # Сначала проверяем, может ID есть в основном результате
            ensembl = hit.get('ensembl')
            if ensembl:
                if isinstance(ensembl, dict):
                    gene_id = ensembl.get('gene')
                    if gene_id:
                        return gene_id
                elif isinstance(ensembl, list) and len(ensembl) > 0:
                    gene_id = ensembl[0].get('gene')
                    if gene_id:
                        return gene_id
            
            # Собираем все синонимы
            aliases = []
            
            # alias может быть строкой или списком
            alias = hit.get('alias', [])
            if isinstance(alias, str):
                aliases.append(alias)
            elif isinstance(alias, list):
                aliases.extend(alias)
            
            # other_names тоже может быть строкой или списком
            other_names = hit.get('other_names', [])
            if isinstance(other_names, str):
                aliases.append(other_names)
            elif isinstance(other_names, list):
                aliases.extend(other_names)
            
            # Убираем дубликаты и исходный символ
            aliases = list(set(aliases))
            aliases = [a for a in aliases if a != gene_symbol]
            
            if not aliases:
                return None
            
            if verbose:
                print(f"  {gene_symbol} → синонимы: {aliases[:5]}...")
            
            # Пробуем искать по каждому синониму
            for syn in aliases:
                syn_query = mg.query(
                    syn,
                    scopes='symbol',
                    fields='ensembl.gene',
                    species='human',
                    size=1
                )
                
                if syn_query and 'hits' in syn_query and syn_query['hits']:
                    ensembl_data = syn_query['hits'][0].get('ensembl')
                    
                    if isinstance(ensembl_data, dict):
                        gene_id = ensembl_data.get('gene')
                        if gene_id:
                            return gene_id
                    elif isinstance(ensembl_data, list) and len(ensembl_data) > 0:
                        gene_id = ensembl_data[0].get('gene')
                        if gene_id:
                            return gene_id
            
            return None
            
        except Exception as e:
            if verbose:
                print(f"  ⚠ Ошибка при поиске синонимов для {gene_symbol}: {e}")
            return None
    
    mg = MyGeneInfo()
    
    # 1) Первичный запрос по списку символов
    query = mg.querymany(
        list_of_genes,
        scopes='symbol',
        fields='ensembl.gene',
        species='human',
        as_dataframe=False
    )
    
    results = {}
    failed_genes = []
    
    for entry in tqdm(query, desc="Primary search"):
        sym = entry.get('query')
        ensembl = entry.get('ensembl')
        notfound = entry.get('notfound', False)
        
        gene_id = None
        
        # Парсим ensembl
        if isinstance(ensembl, dict):
            gene_id = ensembl.get('gene')
        elif isinstance(ensembl, list) and len(ensembl) > 0:
            gene_id = ensembl[0].get('gene')
        
        results[sym] = gene_id
        
        # Если не нашли, добавляем в список для поиска по синонимам
        if gene_id is None or notfound:
            failed_genes.append(sym)
    
    # 2) Поиск по синонимам для ненайденных генов
    if use_synonyms and failed_genes:
        if verbose:
            print(f"\n🔍 Поиск по синонимам для {len(failed_genes)} генов...")
        
        for gene in tqdm(failed_genes, desc="Synonym search"):
            synonym_result = _search_by_synonyms(gene, mg, verbose)
            if synonym_result:
                results[gene] = synonym_result
                if verbose:
                    print(f"  ✓ {gene} → {synonym_result}")
    
    # 3) Статистика
    found = sum(1 for v in results.values() if v is not None)
    total = len(list_of_genes)
    
    print(f"\n📊 Результаты:")
    print(f"  Найдено: {found}/{total} ({100*found/total:.1f}%)")
    print(f"  Не найдено: {total-found}")
    
    return results
