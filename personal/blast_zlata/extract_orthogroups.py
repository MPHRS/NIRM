from Bio import SeqIO
from collections import defaultdict
import os

# Путь к папке с протеомами
proteome_dir = "proteomes"
orthogroup_file = "orthogroups.txt"
output_dir = "orthogroup_fastas"

# Чтение всех белков из протеомов
print("📦 Чтение протеомов...")
all_seqs = {}
for filename in os.listdir(proteome_dir):
    if filename.endswith(".faa"):
        for record in SeqIO.parse(os.path.join(proteome_dir, filename), "fasta"):
            all_seqs[record.id] = record

# Создание выходной папки
os.makedirs(output_dir, exist_ok=True)

print("✂️ Извлечение ортогрупп...")
with open(orthogroup_file) as f:
    for i, line in enumerate(f):
        ids = line.strip().split()
        sequences = [all_seqs[seq_id] for seq_id in ids if seq_id in all_seqs]
        if len(sequences) >= 3:  # минимум 3 вида для дерева
            out_file = os.path.join(output_dir, f"group_{i+1}.faa")
            SeqIO.write(sequences, out_file, "fasta")

print("✅ Готово: FASTA-файлы с ортогруппами лежат в 'orthogroup_fastas/'")

