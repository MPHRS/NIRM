
# align_and_tree.py
#!/usr/bin/env python3
"""
align_and_tree.py: Скрипт для конкатенации CDS porB+tbpB,
выравнивания всех штаммов и построения филогенетического дерева.
Requirements: MAFFT, IQ-TREE, Biopython
Usage:
  python3 align_and_tree.py <input_dir> <concat_fasta> <alignment_fasta> <tree_prefix>
"""
import os
import sys
from Bio import SeqIO
def concat_genes(input_dir, concat_fasta):
    os.makedirs(os.path.dirname(concat_fasta), exist_ok=True)

    records = []
    for fn in os.listdir(input_dir):
        if fn.endswith('.fasta'):
            strain = os.path.splitext(fn)[0]
            path = os.path.join(input_dir, fn)
            seqs = list(SeqIO.parse(path, 'fasta'))
            if len(seqs) != 2:
                print(f"[WARN] {strain}: expected 2 sequences, found {len(seqs)}")
                continue
            porb = [s for s in seqs if 'porb' in s.description.lower()][0]
            tbpb = [s for s in seqs if 'tbpb' in s.description.lower() or 'transferrin-binding' in s.description.lower()][0]
            porb.seq += tbpb.seq
            porb.id = strain
            porb.description = ''
            records.append(porb)

    SeqIO.write(records, concat_fasta, 'fasta')
    print(f"[INFO] Concatenated {len(records)} strains to {concat_fasta}")

def run_cmd(cmd):
    print(f"[RUN] {cmd}")
    r = os.system(cmd)
    if r!=0:
        sys.exit(f"Command failed: {cmd}")

def main():
    if len(sys.argv)!=5:
        print(__doc__)
        sys.exit(1)
    inp, concat_fasta, aln_fasta, tree_prefix = sys.argv[1:]
    concat_genes(inp, concat_fasta)
    run_cmd(f"mafft --auto {concat_fasta} > {aln_fasta}")
    run_cmd(f"iqtree -s {aln_fasta} -nt AUTO -m MFP -bb 1000 -pre {tree_prefix}")
    print(f"[OK] Alignment: {aln_fasta}, Tree files prefix: {tree_prefix}")

if __name__=='__main__':
    main()