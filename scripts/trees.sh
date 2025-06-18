#!/bin/bash

# Параметр: максимальное число топовых генов, для которых строится дерево
max_genes=10

# Создаем каталог для итоговых деревьев, если он не существует
mkdir -p final_trees

# Обходим варианты от 1 до max_genes
for (( N=1; N<=max_genes; N++ )); do
  echo "Обработка варианта с топовыми $N генами"

  variant_dir="trees_temp/trees_${N}_genes_mid"
  mkdir -p "$variant_dir"

  # Путь к готовому конкатенированному файлу
  concat_fasta="tmp_files/median_shenn/variant_${N}/concat.fasta"

  if [ ! -f "$concat_fasta" ]; then
    echo "Ошибка: файл $concat_fasta не найден!" >&2
    continue
  fi

  # Построение дерева с IQ-TREE (автоматический подбор модели + 1000 bootstrap)
  iqtree -s "$concat_fasta" -m MFP -bb 1000 -nt AUTO -pre "${variant_dir}/tree_${N}genes"

  tree_file="${variant_dir}/tree_${N}genes.treefile"
  if [ -f "$tree_file" ]; then
    echo "Дерево для $N генов успешно построено."
    cp "$tree_file" "mid_shen_tree/tree_${N}genes.treefile"
  else
    echo "Ошибка: дерево не построено для $N генов!" >&2
  fi

  echo "------------------------------------------------------------"
done

echo "Все варианты обработаны. Итоговые деревья находятся в каталоге final_trees/"
