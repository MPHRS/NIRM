#!/usr/bin/env python3
import csv
import subprocess
from pathlib import Path
from Bio import SeqIO
import numpy as np
import dendropy

# --- Константы ---
WORKDIR = Path("mlst_extracted")
MAPPING   = Path("tmp_files/mapping.csv")
REF_TREE  = Path("final_trees/tree_all_genes.treefile")

CONCAT_F = WORKDIR/"concat_mlst.fasta"
ALIGNED_F = WORKDIR/"aligned_mlst.fasta"
IQTREE_PREFIX = WORKDIR/"mlst_iqtree"
RESULTS = WORKDIR/"results.txt"

# 1) Читаем mapping.csv
mapping = {}
with open(MAPPING) as f:
    reader = csv.DictReader(f)
    for row in reader:
        acc = row["assembly_accession"]
        # из "strain=WHO_Z_2024" → WHO_Z_2024
        strain = row["strain"].split("=",1)[1]
        mapping[acc] = f"strain_{strain}"

# 2) Список генов и загрузка аллелей
gene_files = sorted(WORKDIR.glob("*.fasta"))
# gene_name → { acc: SeqRecord }
gene_alleles = {}
max_len = {}
for gf in gene_files:
    gene = gf.stem
    recs = list(SeqIO.parse(gf, "fasta"))
    d = {}
    for r in recs:
        # r.id = "GCF_xxx|lcl|...." → берем префикс до '|'
        acc = r.id.split("|",1)[0]
        d[acc] = r.seq
    gene_alleles[gene] = d
    # максимальная длина этого гена
    max_len[gene] = max(len(s) for s in d.values())

# все штаммы по mapping
all_acc = sorted(mapping.keys())

# 3) Конкатенация
with open(CONCAT_F, "w") as outf:
    for acc in all_acc:
        seqs = []
        for gene in gene_files:
            g = gene.stem
            allele = gene_alleles[g].get(acc)
            if allele:
                seq = str(allele)
            else:
                # заполняем gaps нужной длины
                seq = "-" * max_len[g]
            seqs.append(seq)
        full = "".join(seqs)
        outf.write(f">{mapping[acc]}\n{full}\n")

# 4) MAFFT
subprocess.run([
    "mafft", "--auto",
    str(CONCAT_F)],
    stdout=open(ALIGNED_F, "w"),
    check=True
)

# 5) IQ-TREE с 100 bootstrap
subprocess.run([
    "iqtree2",
    "-s", str(ALIGNED_F),
    "-m", "GTR+G",
    "-B", "100",
    "-pre", str(IQTREE_PREFIX)
], check=True)

# 6) RF-distance
ref = dendropy.Tree.get(path=str(REF_TREE), schema="newick")
mlst = dendropy.Tree.get(path=str(IQTREE_PREFIX+".treefile"), schema="newick")
rf = ref.symmetric_difference(mlst)

# 7) Средняя Bootstrap-поддержка
# IQ-TREE в .treefile хранит в узлах числа в формате "{support}:branch_length"
# DendroPy парсит их как node.label
bs = []
for node in mlst.internal_nodes():
    if node.label is not None:
        try:
            bs.append(float(node.label))
        except ValueError:
            pass
mean_bs = np.mean(bs) if bs else 0

# 8) Результаты
with open(RESULTS, "w") as f:
    f.write(f"RF distance: {rf}\n")
    f.write(f"Mean bootstrap support: {mean_bs:.1f}%\n")

print("Done. Results in", RESULTS)
