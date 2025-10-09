import pandas as pd
import gzip

# Если нужны ENS названия 

gtf_path = "path_to/Homo_sapiens.GRCh38.115.gtf.gz"
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

import gzip
import pandas as pd

# Если нужны hgnc названия 
records = []

with gzip.open(gtf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        if feature != flag:
            continue

        try:
            start, end = int(parts[3]), int(parts[4])
        except ValueError:
            continue

        attr_field = parts[8]
        # Устойчивый парсинг атрибутов GTF: key "value";
        attrs = {}
        for item in attr_field.split(";"):
            item = item.strip()
            if not item:
                continue
            kv = item.split(" ", 1)
            if len(kv) != 2:
                continue
            key = kv[0]
            val = kv[1].strip().strip('"')
            attrs[key] = val

        gene_name = attrs.get("gene_name")
        # если по каким‑то причинам gene_name отсутствует, можно падать обратно на gene_id
        if gene_name is None:
            gene_name = attrs.get("gene_id")
        if gene_name:
            records.append((gene_name, start, end))
