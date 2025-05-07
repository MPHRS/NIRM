import os
import subprocess
import pandas as pd
from Bio import SeqIO
import numpy as np
from scipy.stats import entropy

# Файл proteinortho
proteinortho_tsv = "cds_ortho.proteinortho.tsv"
core_dir = "./cds/"
tmp_alignment_file = "tmp_alignment.faa"
mafft_logs_dir = "./mafft_logs/"

# Создаём папку для логов
os.makedirs(mafft_logs_dir, exist_ok=True)

# Парсим proteinortho tsv
df = pd.read_csv(proteinortho_tsv, sep="\t")

# Выбираем только ортогруппы без паралогов и с Alg.-Conn. = 1
#filtered_df = df[(df["Alg.-Conn."] == 1) & (df["Genes"] == len(df.columns[3:]))]
#отбираются ортологичные группы, где количество видов (# Species) и количество генов (Genes) совпадает с количеством используемых образцов, а параметр "Alg.-Conn." не учитывается.
filtered_df = df[(df["# Species"] == len(df.columns[3:])) & (df["Genes"] == len(df.columns[3:]))]

# Печать сколько ортогрупп получилось
print(f"Ортогрупп для анализа: {len(filtered_df)}")

# Словарь для хранения значений Шеннона
shannon_indices = {}

# Проходим по каждой ортогруппе
for idx, row in filtered_df.iterrows():
    sequences = []
    
    # Перебираем образцы (после третьей колонки)
    for sample, genes in row.iloc[3:].items():
        if pd.isna(genes) or genes == "*":
            continue

        gene_id = genes.split(",")[0]  # берем первый (он единственный) ген

        core_file = os.path.join(core_dir, sample + ".core")
        
        # Ищем нужный ген в core-файле
        found = False
        for record in SeqIO.parse(core_file, "fasta"):
            if record.id.split()[0] == gene_id:
                sequences.append(record)
                found = True
                break

        if not found:
            print(f"Не найден ген {gene_id} в файле {core_file}")

    # Записываем последовательности во временный файл
    SeqIO.write(sequences, tmp_alignment_file, "fasta")

    # Выравниваем последовательности при помощи mafft
    aligned_file = f"{idx}_aligned.fna"
    log_file = os.path.join(mafft_logs_dir, f"{idx}_mafft.log")

    # Запускаем MAFFT и сохраняем логи
    with open(aligned_file, "w") as aligned, open(log_file, "w") as log:
        subprocess.run(["mafft", "--auto", tmp_alignment_file], stdout=aligned, stderr=log)

    # Читаем выравнивание 
    alignment = list(SeqIO.parse(aligned_file, "fasta"))
    aln_array = np.array([list(str(rec.seq)) for rec in alignment])

    # Подсчет частот аминокислот по столбцам
    shannon_values = []
    for col in aln_array.T:
        _, counts = np.unique(col, return_counts=True)
        probs = counts / counts.sum()
        shannon_values.append(entropy(probs, base=2))

    # Средний индекс Шеннона по всей ортогруппе
    mean_shannon = np.mean(shannon_values)
    shannon_indices[idx] = mean_shannon

    print(f"Ортогруппа {idx}: Индекс Шеннона = {mean_shannon:.4f}")

# Записываем индексы Шеннона в файл
with open("shannon_indices_cds.tsv", "w") as out:
    out.write("Orthogroup\tShannonIndex\n")
    for idx, val in shannon_indices.items():
        out.write(f"{idx}\t{val:.4f}\n")

# Удаляем временный файл
os.remove(tmp_alignment_file) 