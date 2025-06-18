#!/bin/bash

# Папка с MLST аллелями
MLST_DIR="./db/pubmlst/neisseria"
# Папка с CDS всех штаммов
GENOME_DIR="./data/genomes/refseq/bacteria"
# Папка, куда сохраняем гены
OUTDIR="mlst_genes"

mkdir -p "$OUTDIR"

# Все гены MLST-схемы
GENES=(abcZ adk aroE fumC gdh pdhC pgm)

for gene in "${GENES[@]}"; do
    echo "=== processing $gene ==="

    query="$MLST_DIR/${gene}.tfa"
    out="$OUTDIR/${gene}.fasta"
    > "$out"  # очищаем выходной файл

    for genome_path in "$GENOME_DIR"/GCF_*; do
        sample=$(basename "$genome_path")

        # Ищем соответствующий .fna-файл
        cds=$(find "$genome_path" -name '*cds_from_genomic.fna' | head -n 1)
        if [[ ! -f "$cds" ]]; then
            echo "  !! CDS file not found for $sample"
            continue
        fi

        echo "  >> sample: $sample"
        echo "     using CDS file: $cds"

        # BLAST: ищем лучший хит
        hit=$(blastn -query "$query" -subject "$cds" \
            -perc_identity 85 -qcov_hsp_perc 80 \
            -max_target_seqs 1 -outfmt "6 sseqid sstart send sseq" 2>/dev/null)

        if [[ -z "$hit" ]]; then
            echo "     no match found"
            continue
        fi

        # Парсим BLAST результат
        sseqid=$(echo "$hit" | awk '{print $1}')
        start=$(echo "$hit" | awk '{print $2}')
        end=$(echo "$hit" | awk '{print $3}')
        seq=$(echo "$hit" | awk '{print $4}')

        # Проверяем направление
        if (( start > end )); then
            seq=$(echo "$seq" | tr "ACGTacgt" "TGCAtgca" | rev)
        fi

        # Пишем в FASTA
        echo ">${sample}_${sseqid}" >> "$out"
        echo "$seq" >> "$out"
    done
done

echo "DONE. Extracted gene sequences saved to $OUTDIR"
