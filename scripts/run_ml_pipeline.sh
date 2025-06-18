#!/usr/bin/env bash
set -euo pipefail

# == Параметры по умолчанию ==
TOP_LIST="tmp_files/top_ml_genes.txt"
ALIGN_DIR="alignments"
CONCAT_DIR="concat_ml"
TREES_DIR="final_trees_ml"
TMP_DIR="trees_ml_temp"
BOOTSTRAP_SUMMARY="bootstrap_summary_ml.txt"

# == Опции командной строки ==
usage() {
  echo "Usage: $0 [-t top_list] [-i align_dir] [-c concat_dir] \
             [-r trees_dir] [-m tmp_dir] [-b bootstrap_summary]"
  exit 1
}

while getopts "t:i:c:r:m:b:" opt; do
  case $opt in
    t) TOP_LIST="$OPTARG"       ;;  # файл со списком топ-N генов
    i) ALIGN_DIR="$OPTARG"      ;;  # директория с *.aln
    c) CONCAT_DIR="$OPTARG"     ;;  # куда складывать конкатенаты
    r) TREES_DIR="$OPTARG"      ;;  # куда складывать финальные деревья
    m) TMP_DIR="$OPTARG"        ;;  # временная директория для iqtree
    b) BOOTSTRAP_SUMMARY="$OPTARG" ;; # файл сводки bootstrap
    *) usage                     ;;  # неизвестная опция
  esac
done

# == Проверка существования входного файла ==
if [[ ! -f "$TOP_LIST" ]]; then
  echo "ERROR: не найден список топ-генов: $TOP_LIST" >&2
  exit 1
fi

# Определяем число генов в списке
MAX_GENES=$(wc -l < "$TOP_LIST")

# Создаём директории
mkdir -p "$CONCAT_DIR" "$TREES_DIR" "$TMP_DIR"

# Загружаем OG в массив
mapfile -t ML_OGS < "$TOP_LIST"

# 1) Конкатенация
for N in $(seq 1 "$MAX_GENES"); do
  variant="${CONCAT_DIR}/variant_${N}"
  mkdir -p "$variant"
  echo "[Concat] Top-${N}"
  for ((i=0; i<N; i++)); do
    og=${ML_OGS[$i]}
    src="${ALIGN_DIR}/${og}.aln"
    if [[ -f "$src" ]]; then
      cp "$src" "$variant/"
    else
      echo "WARNING: $src не найден" >&2
    fi
  done

  # Python конкатенация
  python3 - <<EOF
import os
from Bio import AlignIO
variant_dir = "$variant"
files = sorted(f for f in os.listdir(variant_dir) if f.endswith('.aln'))
if not files:
    raise SystemExit(f"No .aln in {variant_dir}")
first = AlignIO.read(os.path.join(variant_dir, files[0]), "fasta")
strains = [r.id for r in first]
concat = {s: "" for s in strains}
for fn in files:
    aln = AlignIO.read(os.path.join(variant_dir, fn), "fasta")
    for rec, s in zip(aln, strains):
        concat[s] += str(rec.seq)
out_fasta = os.path.join(variant_dir, "concat.fasta")
with open(out_fasta, "w") as out:
    for s, seq in concat.items():
        out.write(f">{s}\n{seq}\n")
EOF

done

# 2) IQ-TREE
for N in $(seq 1 "$MAX_GENES"); do
  echo "[IQ-TREE] Top-${N}"
  work="${TMP_DIR}/trees_${N}"
  mkdir -p "$work"
  concat_fa="${CONCAT_DIR}/variant_${N}/concat.fasta"
  if [[ ! -f "$concat_fa" ]]; then
    echo "WARNING: $concat_fa не найден, пропускаем" >&2
    continue
  fi
  iqtree -s "$concat_fa" -m MFP -bb 1000 -nt AUTO -pre "$work/tree" >/dev/null
  cp "$work/tree.treefile" "${TREES_DIR}/tree_${N}.treefile"
done

# 3) Сводка bootstrap
echo -e "Genes\tAverage_Bootstrap" > "$BOOTSTRAP_SUMMARY"
for tf in "${TREES_DIR}"/tree_*.treefile; do
  n=$(basename "$tf" | sed -E 's/tree_([0-9]+)\.treefile/\1/')
  boots=$(grep -oP '\)\K[0-9]+' "$tf" || echo "")
  if [[ -n "$boots" ]]; then
    avg=$(echo "$boots" | awk '{sum+=$1} END{printf "%.2f", sum/NR}')
  else
    avg=0
  fi
  echo -e "${n}\t${avg}" >> "$BOOTSTRAP_SUMMARY"
done

echo "Done: trees in $TREES_DIR, summary in $BOOTSTRAP_SUMMARY"
