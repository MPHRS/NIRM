for file in blast_results/*.txt; do
    awk '!seen[$1]++' "$file" > "${file%.txt}_best.txt"
done

