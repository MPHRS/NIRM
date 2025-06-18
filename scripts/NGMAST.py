#!/usr/bin/env python3
"""
NGMAST.py: Скрипт для извлечения CDS-последовательностей porB и tbpB из оригинальных FASTA-файлов (_cds_from_genomic.fna),
с сохранением заголовка и добавлением strain в заголовок.
Usage:
  python3 NGMAST.py <root_dir> <output_dir>
"""

import os
import sys
import re
from Bio import SeqIO


def extract_strain_name(filename):
    # Извлекает WHO_*_YYYY часть из имени файла
    match = re.search(r'_(WHO_[A-Za-z0-9_]+)_cds_from_genomic\.fna$', filename)
    return match.group(1) if match else None


def main(root_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"[INFO] Created output directory: {output_dir}")

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if not filename.endswith('_cds_from_genomic.fna'):
                continue

            filepath = os.path.join(dirpath, filename)
            strain = extract_strain_name(filename)

            if not strain:
                print(f"[WARNING] Skipping file without valid strain name: {filename}")
                continue

            gene_dict = {}
            for record in SeqIO.parse(filepath, "fasta"):
                header = record.description.lower()
                if 'porb' in header:
                    record.id = f"{record.id}_{strain}"
                    record.description += f" [strain={strain}]"
                    gene_dict['porB'] = record
                # Ищем tbpB: прямое упоминание, полное название, или "protein-like" без tbpA
                elif ('tbpb' in header
                      or ('transferrin-binding protein' in header and 'tbpa' not in header)):
                    record.id = f"{record.id}_{strain}"
                    record.description += f" [strain={strain}]"
                    gene_dict['tbpB'] = record

            if 'porB' in gene_dict and 'tbpB' in gene_dict:
                out_file = os.path.join(output_dir, f"{strain}.fasta")
                SeqIO.write([gene_dict['porB'], gene_dict['tbpB']], out_file, "fasta")
                print(f"[OK] Wrote {out_file}")
            else:
                missing = {'porB', 'tbpB'} - set(gene_dict.keys())
                print(f"[INFO] {filename}: missing genes {missing}, no output generated")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
