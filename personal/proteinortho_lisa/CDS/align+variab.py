import os
import subprocess
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

proteinortho_tsv = "cds_ortho.proteinortho.tsv"
core_dir = "./cds/"
tmp_alignment_file = "tmp_alignment.faa"


df = pd.read_csv(proteinortho_tsv, sep="\t")

# Фильтрация ортогрупп
filtered_df = df[(df["# Species"] == len(df.columns[3:])) & (df["Genes"] == len(df.columns[3:]))]

# Печать сколько ортогрупп получилось
print(f"Ортогрупп для анализа: {len(filtered_df)}")

identity_results = {}

for idx, row in filtered_df.iterrows():
    sequences = []
    for sample, genes in row.iloc[3:].items():
        if pd.isna(genes) or genes == "*":
            continue

        gene_id = genes.split(",")[0]
        core_file = os.path.join(core_dir, sample + ".core")

        found = False
        for record in SeqIO.parse(core_file, "fasta"):
            if record.id.split()[0] == gene_id:
                sequences.append(record)
                found = True
                break

        if not found:
            print(f"Не найден ген {gene_id} в файле {core_file}")

    SeqIO.write(sequences, tmp_alignment_file, "fasta")

    aligned_file = f"{idx}_aligned.fna"
    

    # Выполняем выравнивание через MAFFT с сохранением логов
    subprocess.run(["mafft", "--auto", tmp_alignment_file], stdout=open(aligned_file, "w"))

    alignment = AlignIO.read(aligned_file, "fasta")
    aln_array = np.array([list(rec.seq) for rec in alignment])

    # Подсчёт процента идентичности
    identical_positions = sum(len(set(col)) == 1 for col in aln_array.T)
    total_positions = aln_array.shape[1]

    identity_percent = (identical_positions / total_positions) * 100
    variability_percent = 100 - identity_percent

    identity_results[idx] = (identity_percent, variability_percent)

    print(f"Ортогруппа {idx}: Идентичность = {identity_percent:.2f}%, Вариабельность = {variability_percent:.2f}%")

# Сохраняем результаты в файл
with open("identity_variability_cds.tsv", "w") as out:
    out.write("Orthogroup\tIdentity(%)\tVariability(%)\n")
    for idx, (identity, variability) in identity_results.items():
        out.write(f"{idx}\t{identity:.2f}\t{variability:.2f}\n")

os.remove(tmp_alignment_file)
