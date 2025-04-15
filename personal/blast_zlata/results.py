import os

# Папка с результатами BLAST
blast_results_dir = "blast_results"

# Создадим список для хранения данных
gene_data = []

# Обработаем все файлы в папке с результатами BLAST
for file_name in os.listdir(blast_results_dir):
    if file_name.endswith(".txt"):  # Только файлы с результатами BLAST
        file_path = os.path.join(blast_results_dir, file_name)
        
        with open(file_path, 'r') as file:
            for line in file:
                # Разбираем строку
                parts = line.strip().split("\t")
                
                # Извлекаем идентификаторы, процент идентичности и e-value
                query = parts[0]
                subject = parts[1]
                identity = float(parts[2])
                e_value = float(parts[10])

                # Сохраняем данные в список
                gene_data.append({
                    'query': query,
                    'subject': subject,
                    'identity': identity,
                    'e_value': e_value
                })

# Отсортируем данные по проценту идентичности и e-value
sorted_by_identity = sorted(gene_data, key=lambda x: x['identity'], reverse=True)  # Отсортировать по убыванию идентичности
sorted_by_e_value = sorted(gene_data, key=lambda x: x['e_value'])  # Отсортировать по возрастанию e-value

# Откроем файл для записи результатов
with open('top_identical_genes.txt', 'w') as f:
    f.write("Топ-10 по проценту идентичности:\n")
    for item in sorted_by_identity[:10]:
        f.write(f"{item['query']} - {item['subject']}: {item['identity']}%\n")

with open('top_evalue_genes.txt', 'w') as f:
    f.write("\nТоп-10 по e-value:\n")
    for item in sorted_by_e_value[:10]:
        f.write(f"{item['query']} - {item['subject']}: e-value = {item['e_value']}\n")

print("Результаты записаны в файлы: top_identical_genes.txt и top_evalue_genes.txt")

