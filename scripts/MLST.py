#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
from Bio import SeqIO

# --- Настройки ---
# Путь к папке с геномами
BACT_DIR = Path.home() / "PISH/NIRM/16_05_2025/data/genomes/refseq/bacteria"
# Путь к папке с аллелями MLST Neisseria
MLST_DIR = Path.home() / "miniconda3/envs/GR_data/db/pubmlst/neisseria"
# Папка для результатов
OUT_DIR = Path.cwd() / "mlst_extracted"

# Параметры BLAST
EVALUE_CUTOFF = 1e-5
THREADS = 4

# --- Функции ---
def make_blast_db(cds_fasta: Path):
    """Создать BLAST-базу из fasta-файла CDS."""
    cmd = [
        "makeblastdb",
        "-in", str(cds_fasta),
        "-dbtype", "nucl",
        "-out", str(cds_fasta.with_suffix(""))
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

def blast_and_get_best(cds_db: Path, query_fasta: Path):
    """Запустить blastn и вернуть id лучшего гита (query слева, subject id справа)."""
    outfmt = "6 qseqid sseqid pident length evalue bitscore"
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", str(cds_db),
        "-evalue", str(EVALUE_CUTOFF),
        "-outfmt", outfmt,
        "-max_target_seqs", "1",
        "-num_threads", str(THREADS)
    ]
    proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout.splitlines():
        qid, sid, pid, length, evalue, bits = line.split("\t")
        return sid  # возвращаем subject id
    return None

def extract_seq(cds_fasta: Path, seq_id: str):
    """Извлечь SeqRecord по идентификатору из fasta."""
    for rec in SeqIO.parse(str(cds_fasta), "fasta"):
        if rec.id == seq_id:
            return rec
    return None

# --- Основной ход ---
def main():
    OUT_DIR.mkdir(exist_ok=True)
    # Список генов (tfa-файлы)
    gene_files = list(MLST_DIR.glob("*.tfa"))

    # Проход по каждому MLST-гену
    for gene_f in gene_files:
        gene_name = gene_f.stem  # e.g. "abcZ"
        out_fasta = OUT_DIR / f"{gene_name}.fasta"
        records = []

        # Для каждого штамма
        for strain_dir in sorted(BACT_DIR.glob("GCF_*")):
            cds_f = next(strain_dir.glob("*_cds_from_genomic.fna"))
            db_prefix = cds_f.with_suffix("")  # путь до БД

            # создаём базу (если ещё нет)
            if not (db_prefix.with_suffix(".nhr").exists()):
                make_blast_db(cds_f)

            # BLAST
            best_hit = blast_and_get_best(db_prefix, gene_f)
            if best_hit:
                rec = extract_seq(cds_f, best_hit)
                if rec:
                    rec.id = f"{strain_dir.name}|{rec.id}"
                    rec.description = ""
                    records.append(rec)

        # Сохраняем мультиифасту
        if records:
            SeqIO.write(records, str(out_fasta), "fasta")
            print(f"→ {gene_name}: найдено {len(records)} последовательностей")

if __name__ == "__main__":
    main()
