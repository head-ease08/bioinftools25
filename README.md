# 🧬 Bioinformatics Toolkit

> OOP-based Python toolkit for working with biological sequences, FASTQ filtering, BLAST parsing, and GenBank analysis.

---

## 📖 Overview

This project provides a structured, object-oriented framework for common bioinformatics tasks:

- 🔬 **Biological sequence classes** — DNA, RNA, and protein sequences with polymorphic methods
- 🧹 **FASTQ filtering** — powered by Biopython (`SeqIO`, `SeqUtils`)
- 💥 **BLAST result parsing** — extract best hits from BLAST output
- 🧭 **GenBank neighbor extraction** — find neighboring genes in `.gbk` files

---

## 🏗️ Architecture

### Class Hierarchy

```
BiologicalSequence (ABC)
├── NucleicAcidSequence
│   ├── DNASequence
│   └── RNASequence
└── AminoAcidSequence
```

### 🧩 `BiologicalSequence` (Abstract Base Class)

The root of the hierarchy. Defines a common interface for all biological sequences:

| Feature | Method |
|---|---|
| Length | `len(seq)` |
| Indexing & slicing | `seq[0]`, `seq[1:5]` |
| Print output | `print(seq)`, `repr(seq)` |
| Alphabet validation | `seq.check_alphabet()` |

### 🧪 `NucleicAcidSequence`

Implements `BiologicalSequence` for nucleic acids. All methods are **polymorphic** — they rely on `_alphabet` and `_complement_map` defined in subclasses, with no `if/else` branching on sequence type.

| Method | Description |
|---|---|
| `complement()` | Returns complement sequence |
| `reverse()` | Returns reversed sequence |
| `reverse_complement()` | Returns reverse complement |
| `check_alphabet()` | Validates nucleotide alphabet |

> ⚠️ `NucleicAcidSequence` is not meant to be instantiated directly — calling methods on a raw instance raises `NotImplementedError`.

### 🧬 `DNASequence`

Inherits from `NucleicAcidSequence`. Alphabet: `A, T, G, C`.

| Method | Description |
|---|---|
| `transcribe()` | Returns an `RNASequence` (T → U) |

### 🧫 `RNASequence`

Inherits from `NucleicAcidSequence`. Alphabet: `A, U, G, C`.

No additional public methods — everything is inherited.

### 🥩 `AminoAcidSequence`

Implements `BiologicalSequence` for protein sequences. Alphabet: 20 standard amino acids.

| Method | Description |
|---|---|
| `check_alphabet()` | Validates amino acid alphabet |
| `hydrophobic_ratio()` | Returns fraction of hydrophobic residues (A, V, L, I, M, F, W, P) |

---

## 🧹 FASTQ Filtering

The `filter_fastq` function uses **Biopython** to filter reads by three criteria:

| Parameter | Default | Description |
|---|---|---|
| `gc_bounds` | `(0, 100)` | Min/max GC content (%) |
| `length_bounds` | `(0, 2³²)` | Min/max read length |
| `quality_threshold` | `0` | Minimum mean Phred33 quality |

```python
filter_fastq("reads.fastq", "output.fastq",
             gc_bounds=(40, 60),
             length_bounds=(50, 300),
             quality_threshold=20)
```

Filtered reads are saved to the `filtered/` directory.

---

## 💥 BLAST Parser

`parse_blast_results(input_file, output_file)` extracts the **best hit** for each query from a BLAST `.txt` output and saves unique protein descriptions to a file.

---

## 🧭 GenBank Neighbor Extraction

`extract_neighbor_genes(input_gbk, genes, output_fasta, n_before=1, n_after=1)` finds CDS features matching the given gene names and extracts protein sequences of their neighboring genes into a FASTA file.

---

## 📦 Dependencies

- Python 3.8+
- [Biopython](https://biopython.org/) (`Bio.SeqIO`, `Bio.SeqUtils`)

```bash
pip install biopython
```

---

## 🚀 Quick Start

```python
from main import DNASequence, RNASequence, AminoAcidSequence

# 🧬 DNA
dna = DNASequence("ATGCGATC")
print(dna.complement())          # TACGCTAG
print(dna.reverse_complement())  # GATCGCAT
print(dna.transcribe())          # AUGCGAUC

# 🧫 RNA
rna = RNASequence("AUGCGAUC")
print(rna.complement())          # UACGCUAG

# 🥩 Protein
protein = AminoAcidSequence("ACDEFGHIKLMNPQRSTVWY")
print(protein.check_alphabet())  # True
print(protein.hydrophobic_ratio())

# 🧹 FASTQ filtering
from main import filter_fastq
filter_fastq("input.fastq", "output.fastq",
             gc_bounds=(40, 60), quality_threshold=20)
```

---

## 📁 Project Structure

```
.
├── main.py              # All classes, functions, and entry points
├── requirements.txt      
├── data/                # Input directory
├── filtered/            # Output directory for filtered FASTQ reads
└── README.md
```

---

## 📝 License

Educational project.