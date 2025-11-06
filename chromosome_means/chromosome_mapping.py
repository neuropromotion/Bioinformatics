import gzip
import re
import pandas as pd

gtf = "annotation.gtf.gz"  # URL for download (GTF): https://www.gencodegenes.org/human/ 
chrs= list(map(str,list(range(1,23)))) # take only 1-22
 
attr_re = re.compile(r'(\S+)\s+"([^"]+)"')

rows = []
with gzip.open(gtf, "rt") as fh:
    for line in fh:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        seqname, source, feature, start, end, score, strand, frame, attributes = parts
        if feature != "gene":
            continue
        attrs = dict(attr_re.findall(attributes))
        #gene_id = attrs.get("gene_id")
        gene_name = attrs.get("gene_name") or attrs.get("Name")
        chrom = seqname[3:]
        if chrom in chrs:
            if gene_name:
                rows.append((gene_name, chrom))

df = pd.DataFrame(rows, columns=["HGNC symbol", "Chromosome/scaffold name"])
df.to_csv("gene_mapping.txt", sep='\t', index=False)
