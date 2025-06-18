#!/usr/bin/env python3
"""
generate_features.py

Проходит по всем MSA (.aln) в указанной папке и для каждого OG вычисляет:
  - Shannon entropy (H)
  - Среднюю длину без gaps
  - SNP-density
  - Число уникальных аллелей
  - Число парсимониально-информативных позиций
  - Среднее попарное расстояние (нормированное)
  - Долю gaps

Сохраняет таблицу (OG × признаки) в CSV.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from Bio import AlignIO


VALID_NT = {"A", "C", "G", "T"}


def compute_stats(aln_path):
    """Вычислить все признаки для одного выравнивания в FASTA (.aln)."""
    aln = AlignIO.read(aln_path, "fasta")
    seqs = [str(rec.seq) for rec in aln]
    nseq = len(seqs)
    L = len(seqs[0])

    # Shannon entropy по колонкам
    ent = []
    for i in range(L):
        col = [s[i].upper() for s in seqs]
        freqs = {b: col.count(b) for b in VALID_NT if col.count(b) > 0}
        total = sum(freqs.values())
        if total == 0:
            ent.append(0.0)
        else:
            p = np.array(list(freqs.values())) / total
            ent.append(-np.sum(p * np.log2(p)))
    H = float(np.mean(ent))

    # Средняя длина без gaps
    lengths = [len(s.replace('-', '')) for s in seqs]
    mean_len = float(np.mean(lengths))

    # SNP-density: доля колонок с энтропией > 0
    snp_density = float(sum(1 for x in ent if x > 0) / L)

    # Количество уникальных аллелей (полные последовательности)
    n_alleles = int(len(set(seqs)))

    # Парсимониально-информативные позиции
    pis = 0
    for i in range(L):
        col = [s[i].upper() for s in seqs]
        counts = [col.count(b) for b in VALID_NT if col.count(b) > 0]
        if sum(1 for cnt in counts if cnt >= 2) >= 2:
            pis += 1
    parsimony_sites = pis

    # Среднее попарное расстояние (нормированная Hamming)
    dists = []
    for i in range(nseq):
        for j in range(i + 1, nseq):
            mismatches = sum(1 for a, b in zip(seqs[i], seqs[j]) if a != b)
            dists.append(mismatches / L)
    mean_pw_dist = float(np.mean(dists)) if dists else 0.0

    # Доля gaps в выравнивании
    total_gaps = sum(s.count('-') for s in seqs)
    gap_fraction = float(total_gaps / (nseq * L))

    return {
        "entropy": H,
        "mean_len": mean_len,
        "snp_density": snp_density,
        "n_alleles": n_alleles,
        "parsimony_sites": parsimony_sites,
        "mean_pairwise_dist": mean_pw_dist,
        "gap_fraction": gap_fraction
    }


def main(args):
    rows = []
    for fn in sorted(os.listdir(args.align_dir)):
        if not fn.endswith('.aln'):
            continue
        og = fn[:-4]  # убрать '.aln'
        path = os.path.join(args.align_dir, fn)
        stats = compute_stats(path)
        stats["OG"] = og
        rows.append(stats)

    if not rows:
        print(f"[ERROR] В папке {args.align_dir} нет файлов .aln", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows).set_index("OG")
    df.to_csv(args.out_csv)
    print(f"✅ Записано {args.out_csv}, shape={df.shape}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Генерация признаков по MSA для каждого OG"
    )
    parser.add_argument(
        "--align-dir", "-i",
        required=True,
        help="Папка с .aln файлами"
    )
    parser.add_argument(
        "--out-csv", "-o",
        default="features.csv",
        help="Путь до выходного CSV"
    )
    args = parser.parse_args()
    main(args)
