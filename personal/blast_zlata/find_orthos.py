from Bio import SeqIO
import os

# Путь к папке с протеомами
proteome_dir = "/home/zlata/grant/ncbi_dataset/data"  # Укажите путь к папке с протеомами

# Путь к файлу с ортогруппами
orthogroup_file = "orthogroups.txt"  # Если файл в той же директории, просто имя

# Папка для сохранения последовательностей ортогрупп
output_dir = "orthogroup_fastas"

# Чтение всех белков из протеомов
print("📦 Чтение протеомов...")
all_seqs = {}
for subdir, dirs, files in os.walk(proteome_dir):  # os.walk для обхода всех подкаталогов
    for filename in files:
        if filename.endswith(".faa"):  # Протеины в формате .faa
            proteome_file_path = os.path.join(subdir, filename)
            for record in SeqIO.parse(proteome_file_path, "fasta"):
                all_seqs[record.id] = record

# Создание выходной папки
os.makedirs(output_dir, exist_ok=True)

print("✂️ Извлечение ортогрупп...")
# Открытие файла с ортогруппами с обработкой ошибок кодировки
with open(orthogroup_file, encoding='ISO-8859-1', errors='ignore') as f:  # Можно использовать 'ISO-8859-1' или 'windows-1251'
    for i, line in enumerate(f):
        ids = line.strip().split()
        sequences = [all_seqs[seq_id] for seq_id in ids if seq_id in all_seqs]
        if len(sequences) >= 3:  # минимум 3 вида для дерева
            out_file = os.path.join(output_dir, f"group_{i+1}.faa")
            SeqIO.write(sequences, out_file, "fasta")

print(f"✅ Готово: FASTA-файлы с ортогруппами сохранены в '{output_dir}'")

