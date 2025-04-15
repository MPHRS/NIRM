from Bio import SeqIO
import os

trim_dir = "/home/zlata/grant/ncbi_dataset/trimmed"  # Путь к папке с триммированными файлами
out_file = "concatenated_alignment.faa"  # Имя итогового файла

genomes = set()
seqs_per_genome = {}

# Обрабатываем каждый файл по очереди
with open(out_file, "w") as out:
    for filename in sorted(os.listdir(trim_dir)):
        if filename.endswith(".trim.faa"):
            file_path = os.path.join(trim_dir, filename)
            records = list(SeqIO.parse(file_path, "fasta"))
            id_map = {r.id.split('.')[0]: r for r in records}  # Геном = до первой точки
            
            for genome in id_map:
                genomes.add(genome)
            
            for g in genomes:
                if g in id_map:
                    # Записываем последовательности в итоговый файл по одному
                    out.write(f">{g}\n{str(id_map[g].seq)}\n")
                else:
                    # Заполнение пропущенных групп gap-ами
                    length = len(records[0].seq)
                    out.write(f">{g}\n{'-' * length}\n")

