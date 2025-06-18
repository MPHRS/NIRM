#!/usr/bin/env python3
from Bio import SeqIO
import csv
import sys
import os

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <cds_fasta_file.fna> <mapping.csv>")
    sys.exit(1)

fasta_file = sys.argv[1]
mapping_file = sys.argv[2]

# Читаем mapping: ключ – assembly accession, значение – (biosample, strain)
mapping = {}
with open(mapping_file, 'r') as mf:
    reader = csv.DictReader(mf)
    for row in reader:
        acc = row["assembly_accession"]
        biosample = row["biosample"]
        strain = row["strain"]
        mapping[acc] = (biosample, strain)

# Извлекаем assembly accession из имени файла
# Пример: GCF_001047225.1_ASM104722v1_cds_from_genomic.fna -> GCF_001047225.1
basename = os.path.basename(fasta_file)
assembly = basename.split("_")[0] + "_" + basename.split("_")[1]

# Генерируем имя выходного файла
output_file = os.path.splitext(fasta_file)[0] + "_annotated.fa"

records = []
for record in SeqIO.parse(fasta_file, "fasta"):
    if assembly in mapping:
        biosample, strain = mapping[assembly]
        # Добавляем к ID аннотацию
        record.id = f"{record.id}|{biosample}|{strain}"
        record.description = ""  # очищаем description
    else:
        record.id = record.id  # оставляем без изменений
    records.append(record)

# Запись аннотированного файла
SeqIO.write(records, output_file, "fasta")
print(f"Annotated CDS saved to {output_file}")