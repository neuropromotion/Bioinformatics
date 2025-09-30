import pandas as pd
import gzip

# Укажите путь к GTF-файлу Ensembl (распакованный или .gz)
gtf_path = "/mnt/jack-5/amismailov/miRNA_study/Homo_sapiens.GRCh38.115.gtf.gz"
flag = 'gene'
records = []
with gzip.open(gtf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        feature = parts[2]
        if feature != flag:
            continue
        start, end = int(parts[3]), int(parts[4])
        attr_field = parts[8]
        # Парсим атрибуты GTF
        attrs = {
            item.split(" ")[0]: item.split(" ")[1].strip('";')
            for item in attr_field.split(";") if item.strip()
        }
        gene_id = attrs.get("gene_id")
        if gene_id:
            records.append((gene_id, start, end))
