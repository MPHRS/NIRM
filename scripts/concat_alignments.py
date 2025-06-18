import os
from Bio import AlignIO, SeqIO, SeqRecord, Seq
import pandas as pd

strain_df = pd.read_csv("./tmp_files/mapping.csv")

strains = list(strain_df["strain"] + " " + strain_df["assembly_accession"])

for path in os.listdir('./tmp_files/median_shenn'):
    full_path = os.path.join('./tmp_files/median_shenn', path)
    if not os.path.isdir(full_path):
        continue

    print(f"Processing group: {path}")

    concat_dict = {strain: "" for strain in strains}

    for aligne in os.listdir(full_path):
        if aligne.endswith('.aln'):
            aln_path = os.path.join(full_path, aligne)
            alignment = AlignIO.read(aln_path, "fasta")

            print(f"  Reading {aligne}")

            if len(alignment) != len(strains):
                print(f"  Warning: {aligne} -> {len(alignment)} seqs, expected {len(strains)}")

            for i, record in enumerate(alignment):
                strain_name = strains[i]
                concat_dict[strain_name] += str(record.seq)

    # Создаем SeqRecord'ы
    concat_seqs = []
    for strain, sequence in concat_dict.items():
        concat_seqs.append(SeqRecord.SeqRecord(Seq.Seq(sequence), id=strain, description=""))

    # Путь сохранения concat.fasta внутри папки группы
    output_path = os.path.join(full_path, "concat.fasta")
    with open(output_path, "w") as output_handle:
        SeqIO.write(concat_seqs, output_handle, "fasta")

    print(f"  Saved {output_path}\n")
