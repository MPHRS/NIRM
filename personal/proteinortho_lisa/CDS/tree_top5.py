import pandas as pd
import subprocess
from Bio import SeqIO
import os

# Читаем файл с отфильтрованными индексами Шеннона
df = pd.read_csv('variable_orthogroups_filt.tsv', sep='\t')#('shannon_indices_filt.tsv', sep='\t')
df_sorted = df.sort_values(by= 'Variability(%)', ascending=False).head(5)#(by='ShannonIndex', ascending=False).head(5)  # Берем топ-5

# Путь к директории с файлами ортогрупп
orthogroups_dir = "./mafft_out_v" #"./mafft_out"
# Путь для сохранения результатов
results_dir = "./trees_5_variab/" #"./trees_5_shannon/"
os.makedirs(results_dir, exist_ok=True)

# Функция объединения последовательностей
def concatenate_sequences(group_ids, output_file):
    concatenated = {}
    for group_id in group_ids:
        file_path = f"{orthogroups_dir}/{group_id}_aligned.fna"
        records = list(SeqIO.parse(file_path, "fasta"))
        for record in records:
            if record.id not in concatenated:
                concatenated[record.id] = str(record.seq)
            else:
                concatenated[record.id] += str(record.seq)

    with open(output_file, "w") as output:
        for seq_id, seq in concatenated.items():
            output.write(f">{seq_id}\n{seq}\n")

# Итеративное добавление только 5 ортогрупп
selected_groups = []
for idx, row in df_sorted.iterrows():
    selected_groups.append(int(row["Orthogroup"]))
    print(f"Processing orthogroups: {selected_groups}")

    # Конкатенация текущих ортогрупп
    concatenated_file = f"{results_dir}/concatenated_{len(selected_groups)}.fna"
    concatenate_sequences(selected_groups, concatenated_file)

    # Выравнивание с MAFFT
    aligned_file = f"{results_dir}/aligned_{len(selected_groups)}.fna"
    subprocess.run(f"mafft --auto {concatenated_file} > {aligned_file}", shell=True)

    # Построение дерева
    subprocess.run(f"cd {results_dir} && iqtree -s aligned_{len(selected_groups)}.fna -m GTR+G -bb 1000 -nt AUTO", shell=True)


print("Analysis complete. Top 5 orthogroups processed.")