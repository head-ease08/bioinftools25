# Bioinformatics Command-Line Toolbox

This project is a command-line tool for performing common bioinformatics tasks, including DNA/RNA sequence manipulation and filtering of FASTQ data.

## 🚀 Features

The script has two main modes of operation:

1.  **Simple Sequence Processing (`process`)**:
    *   Accepts one or more sequences as string inputs.
    *   Performs one of the following operations:
        *   `transcribe`: Transcribes DNA into RNA (T → U).
        *   `reverse`: Reverses the sequence.
        *   `complement`: Creates the complementary strand.
        *   `reverse_complement`: Creates the reverse-complementary strand.
        *   `is_palindrome`: Checks if a sequence is a palindrome.

2.  **FASTQ File Filtering (`filter`)**:
    *   Reads sequences from a FASTQ format file.
    *   Filters them based on specified criteria:
        *   **GC Content**: Keeps sequences within a specified GC content range (in percent).
        *   **Length**: Keeps sequences within a specified length range.
        *   **Quality**: Keeps sequences with an average Phred33 quality score above or equal to a given threshold.
    *   Saves the filtered sequences in FASTA format.

## ⚙️ Installation and Project Structure

This script requires **Python 3.6+**. No external libraries need to be installed. Don't forget to make `main.py` executable!

Your project should have the following structure for the imports from the `modules` directory to work correctly:

├── biotools.py # The main executable script
├── modules/
│ ├── init.py # (Can be empty)
│ ├── is_nucleic_acid.py
│ ├── transcribe.py
│ ├── reverse.py
│ ├── raw_toolbox.py
│ ├── complement_and_reverse.py
│ ├── is_palindrome.py
│ ├── calculate_quality.py
│ └── calculate_gc_content.py
└── README.md # This file

All commands are run from the terminal. The basic syntax is:

```bash
python main.py [arguments]
```

## Simple Sequence Processing

This project contains toolbox for quick operations on sequence strings:

```bash
python main.py [*sequences] [tool]
```

| Tool                 | Description                                  |
| -------------------- | -------------------------------------------- |
| `transcribe`         | Transcribes DNA into RNA (replaces T with U) |
| `reverse`            | Reverses the sequence                        |
| `complement`         | Creates the complementary strand             |
| `reverse_complement` | Creates the reverse-complementary strand     |
| `is_palindrome`      | Checks if the sequence is a palindrome       |
Examples:

# Transcribe a single sequence
python main.py GCTAGTCA transcribe

# Get the reverse-complementary strand
python main.py CGTAGTCAGTCGTATGCGTGTATGCATGTGCATTCATCGATGCATTATTACTATCGGA reverse_complement 

# Check multiple sequences for being a palindrome
python main.py GATATC GAATTC is_palindrome 

FASTQ file filtering: coming soon!

| Argument              | Description                                                                                                                              | Usage Example               |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | --------------------------- |
| `input_file`          | (Required) Path to the input FASTQ file.                                                                                                 | `my_reads.fastq`            |
| `--output_file`       | (Optional) Path to the output FASTA file. If not provided, the result is printed to the console.                                        | `--output_file good.fasta`  |
| `--gc_bounds`         | (Optional) GC content range in percent (min max). Default: `0 100`. A single value can be passed to set only the upper bound.           | `--gc_bounds 40 60`         |
| `--length_bounds`     | (Optional) Sequence length range (min max). Default: `0 4294967296`. A single value can be passed to set only the upper bound.        | `--length_bounds 100 150`   |
| `--quality_threshold` | (Optional) Minimum average Phred33 quality score. Default: `0`.                                                                          | `--quality_threshold 25`    |