# parse_metadata.py
import csv
import sys
import os

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <assembly_summary.txt> <genomes_root_dir>")
    sys.exit(1)

summary_path = sys.argv[1]
genomes_root = sys.argv[2]
output_path = "./tmp_files/mapping.csv"

# 1. Собираем список GCF_... из названий файлов
my_assemblies = set()
for entry in os.listdir(genomes_root):
    if entry.startswith("GCF_") and os.path.isdir(os.path.join(genomes_root, entry)):
        my_assemblies.add(entry)


print(f"Found {len(my_assemblies)} assemblies in genomes/")

# 2. Читаем assembly_summary и фильтруем только нужное
with open(summary_path) as fin, open(output_path, "w", newline="") as fout:
    writer = csv.writer(fout)
    writer.writerow(["assembly_accession", "biosample", "strain"])

    for line in fin:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 20:
            continue
        accession = parts[0]
        biosample = parts[5]
        strain = parts[8]

        # Проверяем, есть ли этот accession (или его вариант) в нашем списке
        if any(accession in a or a.startswith(accession) for a in my_assemblies):
            writer.writerow([accession, biosample, strain])

