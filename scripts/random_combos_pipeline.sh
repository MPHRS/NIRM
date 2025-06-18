#!/usr/bin/env bash
# random_combos_pipeline.sh
# Запускайте из ~/PISH/NIRM/16_05_2025/scripts

set -euo pipefail

# Пути
SCRIPTS_DIR="$(pwd)"
BASE_DIR="${HOME}/PISH/NIRM/16_05_2025"
TMP="${BASE_DIR}/tmp_files"
ALIGN_DIR="${BASE_DIR}/alignments"
MAPPING="${TMP}/mapping.csv"

# Параметры
TOP_N=200      # сколько самых вариабельных брать
K=4            # размер каждой комбинации
REPS=3         # сколько случайных комбинаций на каждый OG
OUT_DIR="${BASE_DIR}/random_combos"
LOG_LIST="${OUT_DIR}/combos_list.txt"

# Подготовка
mkdir -p "${OUT_DIR}"
> "${LOG_LIST}"

echo "=== Начало генерации случайных комбинаций ==="
echo "Топ-${TOP_N} из ${TMP}/sorted_entropy.txt, размер комбо=${K}, по ${REPS} повт."

# Получаем список топ-N OG без расширения .aln
mapfile -t TOP_OGS < <(
  head -n${TOP_N} "${TMP}/sorted_entropy.txt" \
    | awk '{print $1}' \
    | xargs -n1 basename \
    | sed 's/\.aln$//'
)

# Для каждой OG генерируем REPS комбинаций вида [OG + 3 случайных других]
for OG in "${TOP_OGS[@]}"; do
  for r in $(seq 1 "${REPS}"); do
    # выбираем 3 других случайно (без текущего OG)
    OTHERS=( $(printf "%s\n" "${TOP_OGS[@]}" | grep -v "^${OG}$" | shuf -n $((K-1))) )
    COMBO=("${OG}" "${OTHERS[@]}")
    DIR_NAME=$(IFS=_; echo "${COMBO[*]}")
    DIR_PATH="${OUT_DIR}/${DIR_NAME}"
    mkdir -p "${DIR_PATH}"

    echo
    echo ">>> Комбинация: ${DIR_NAME}"
    echo "    Ортогруппы: ${COMBO[*]}"
    echo "    Папка: ${DIR_PATH}"

    # логируем
    echo -e "${DIR_NAME}\t${COMBO[*]}" >> "${LOG_LIST}"

    # копируем .aln файлы
    echo "    Шаг 1/3: копирование .aln..."
    for g in "${COMBO[@]}"; do
      cp "${ALIGN_DIR}/${g}.aln" "${DIR_PATH}/"
      echo "      - ${g}.aln"
    done

    # конкатенация в concat.fasta
    echo "    Шаг 2/3: конкатенация выравниваний в concat.fasta..."
    python3 - <<EOF
import os, csv
from Bio import AlignIO, SeqIO, SeqRecord, Seq

d = "${DIR_PATH}"
# читаем список штаммов
strain_list = []
with open("${MAPPING}") as mf:
    reader = csv.DictReader(mf)
    for row in reader:
        strain_list.append(f"{row['strain']} {row['assembly_accession']}")

# инициализация
concat = {s: "" for s in strain_list}

for fn in sorted(os.listdir(d)):
    if not fn.endswith(".aln"): continue
    path = os.path.join(d, fn)
    aln = AlignIO.read(path, "fasta")
    for rec, strain in zip(aln, strain_list):
        concat[strain] += str(rec.seq)

# запись
out = os.path.join(d, "concat.fasta")
records = [SeqRecord.SeqRecord(Seq.Seq(seq), id=strain, description="")
           for strain, seq in concat.items()]
SeqIO.write(records, out, "fasta")
print(f"Конкатенация готова: {out}")
EOF

    # строим дерево с IQ-TREE
    echo "    Шаг 3/3: запуск IQ-TREE..."
    (
      cd "${DIR_PATH}"
      iqtree -s concat.fasta -m MFP -bb 1000 -nt AUTO -pre tree > tree.log 2>&1
    )
    echo "      IQ-TREE завершён, лог в tree.log"

  done
done

echo
echo "=== Генерация комбинаций завершена ==="
echo "Список в ${LOG_LIST}"
