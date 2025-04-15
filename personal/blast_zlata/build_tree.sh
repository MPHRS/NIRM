#!/bin/bash

set -e

INPUT_DIR="orthogroup_fastas"
ALIGN_DIR="aligned"
TRIM_DIR="trimmed"
mkdir -p "$ALIGN_DIR" "$TRIM_DIR"

echo "🔄 Выравнивание MAFFT..."
for file in $INPUT_DIR/*.faa; do
    base=$(basename "$file" .faa)
    mafft --auto "$file" > "$ALIGN_DIR/$base.aln.faa"
done

echo "✂️ Тримминг выравниваний (trimAl)..."
for file in $ALIGN_DIR/*.aln.faa; do
    base=$(basename "$file" .aln.faa)
    trimal -in "$file" -out "$TRIM_DIR/$base.trim.faa" -automated1
done

echo "🔗 Конкатенация..."
python3 <<EOF
from Bio import SeqIO
import os
from collections import defaultdict

trim_dir = "$TRIM_DIR"
out_file = "concatenated_alignment.faa"

genomes = set()
seqs_per_genome = defaultdict(list)

# Чтение всех выравниваний и добавление по геномам
for filename in sorted(os.listdir(trim_dir)):
    if filename.endswith(".trim.faa"):
        records = list(SeqIO.parse(os.path.join(trim_dir, filename), "fasta"))
        id_map = {r.id.split('.')[0]: r for r in records}  # Genome ID = до первой точки
        for genome in id_map:
            genomes.add(genome)
        for g in genomes:
            if g in id_map:
                seqs_per_genome[g].append(str(id_map[g].seq))
            else:
                # Заполнение пропущенных групп gap-ами
                length = len(records[0].seq)
                seqs_per_genome[g].append('-' * length)

# Запись итогового файла
with open(out_file, "w") as out:
    for genome, chunks in seqs_per_genome.items():
        out.write(f">{genome}\n{''.join(chunks)}\n")
EOF

echo "🌳 Построение дерева (FastTree)..."
FastTree -lg concatenated_alignment.faa > tree.nwk

echo "🎉 Готово! Дерево сохранено в tree.nwk"

