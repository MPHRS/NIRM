#!/bin/bash

# Создаем папку concat_align, если её нет
mkdir -p concat_align

# Извлекаем список .aln файлов из sorted_entropy.txt в порядке убывания вариабельности
readarray -t top_groups < <(awk '{print $1}' ./tmp_files/median_shenn/median_list.txt)

# Создаем подпапки для каждого варианта (1-9 генов) и копируем выравнивания
for N in {1..10}; do
    variant_dir="tmp_files/median_shenn/variant_${N}"
    mkdir -p "$variant_dir"
    
    # Копируем первые N выравниваний
    for ((i=0; i<N; i++)); do
        aln_file="${top_groups[$i]}"
        if [[ -f "$aln_file" ]]; then
            cp "$aln_file" "$variant_dir/"
        else
            echo "Warning: File $aln_file not found!" >&2
        fi
    done
done

echo "Top alignments have been copied to concat_align/variant_{1..10}"