import os

def load_best_hits(filename):
    best_hits = {}
    with open(filename) as f:
        for line in f:
            cols = line.strip().split("\t")
            query, subject = cols[0], cols[1]
            if query not in best_hits:
                best_hits[query] = subject
    return best_hits

blast_dir = "blast_results"
files = sorted([f for f in os.listdir(blast_dir) if f.endswith(".txt")])

rbh_pairs = set()

for file in files:
    file_name = file.replace(".txt", "").replace("_best", "")  # Убираем _best
    parts = file_name.split("_GCA_")

    if len(parts) != 2:
        print(f"⚠️ Пропущен файл {file} (неправильное имя)")
        continue

    org1, org2 = parts[0], "GCA_" + parts[1]

    hits1 = load_best_hits(os.path.join(blast_dir, file))
    reverse_file = f"{org2}_{org1}_best.txt"
    
    if reverse_file not in files:
        print(f"⚠️ Пропущен файл {reverse_file} (нет пары)")
        continue

    hits2 = load_best_hits(os.path.join(blast_dir, reverse_file))

    for query, subject in hits1.items():
        if subject in hits2 and hits2[subject] == query:
            rbh_pairs.add((org1, query, org2, subject))

# Сохранение результатов
with open("reciprocal_best_hits.txt", "w") as out:
    for org1, query, org2, subject in sorted(rbh_pairs):
        out.write(f"{org1}\t{query}\t{org2}\t{subject}\n")

print("✅ Найдены взаимные ортологи! Результаты в reciprocal_best_hits.txt")

