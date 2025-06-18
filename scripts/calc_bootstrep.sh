#!/bin/bash

output_file="bootstrap_summary.txt"
echo -e "Genes\tAverage_Bootstrap" > $output_file

for treefile in final_trees/tree_*genes.treefile; do
    # Вытаскиваем количество генов из имени файла
    genes=$(echo "$treefile" | grep -oP 'tree_\K[0-9]+(?=genes)')

    # Получаем все bootstrap значения (предполагаем формат поддержки: [100] или 100 внутри скобок)
    boots=$(grep -oP '\)[0-9]+' "$treefile" | tr -d ')' )

    # Считаем среднее
    if [ -n "$boots" ]; then
        avg=$(echo "$boots" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print 0}')
    else
        avg=0
    fi

    echo -e "${genes}\t${avg}" >> $output_file
done