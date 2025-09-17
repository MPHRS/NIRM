# NIRMА — Neisseria gonorrhoeae ML-driven Phylogenetics

This repository contains code and analysis from a 6-month research project (grant-supported, ITMO University) on **strain-level phylogenetics and antibiotic resistance profiling of *Neisseria gonorrhoeae***.

---

## Overview

- Built a **pipeline** combining:
  - ortholog detection (OrthoFinder, Proteinortho, BLAST),
  - phylogenetic reconstruction (IQ-TREE2),
  - **machine learning**–based locus selection (Random Forest, DNABERT embeddings),
  - resistance gene screening (CARD database).
- **Input**: FASTA genomes (NCBI RefSeq, WHO-2024 reference set).  
- **Output**: reduced marker panels for high-resolution phylogenetics + ARG profiles.  
- Goal: **resource-efficient, accurate surveillance** of *N. gonorrhoeae* with direct clinical/epidemiological relevance.

 For full methods & results, see [Neisseria gonorrhoeae.pdf](<Neisseria gonorrhoeae.pdf>).  

---

## Repository structure

├── personal/ # individual contributor work (exploratory)   
│ ├── blast_zlata/ # BLAST + CARD outputs     
│ ├── proteinortho_lisa/ # Proteinortho-based analysis       
│ └── seraph/ # OrthoFinder notebooks            
├── scripts/ # final reusable scripts            
├── reports/ # usage instructions (see main.md)              
└── README.md                 

---

## Methods & Tools

- **Input**: FASTA assemblies (NCBI RefSeq, WHO-2024 reference set).  
- **Orthology & alignment**: OrthoFinder, Proteinortho, MAFFT.  
- **Phylogenetics**: IQ-TREE2, RF distance, bootstrap stability.  
- **Machine learning**:  
  - Random Forest on locus features,  
  - DNABERT embeddings for orthogroup pairs.  
- **Resistance screening**: CARD / RGI (`perfect` + `strict` modes).  
- **Languages & libraries**: Python (Biopython, pandas, seaborn, matplotlib), Bash.


---

## Contributors

- Serafim Dobrovolskii (@MPHRS)  
- Kochubei Lisa (@ktsel)  
- Zlata (@Zlata-LAB21)  

---


## Notes

- This repository is primarily maintained as a collaborative space for contributors and as part of the project portfolio.  
- Scripts and pipeline components are not packaged as a ready-to-run tool; see [`reports/main.md`](./reports/main.md) for execution details.

---
## License

[MIT License](./LICENSE)
