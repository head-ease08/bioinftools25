import sys
import os
import argparse
import logging
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

logger = logging.getLogger(__name__)

def parse_blast_results(input_file: str, output_file: str) -> str:
    """
    Reads BLAST txt results and extracts the best hit (first in alignments block)
    for each query. Saves unique protein descriptions in one-column file.
    """
    best_hits = []

    with open(input_file, "r") as infile:
        lines = infile.readlines()

    for i, line in enumerate(lines):
        if line.strip().startswith("Sequences producing significant alignments:"):
            if i + 1 < len(lines):
                first_hit = lines[i + 1].strip()
                if first_hit:
                    parts = first_hit.split(maxsplit=1)
                    if len(parts) > 1:
                        description = parts[1]
                    else:
                        description = parts[0]
                    best_hits.append(description)
    # LOGGING  
    if not best_hits:
        logger.error("No best hits found in %s — output will be empty", input_file)                 
    logger.info("Found %d unique best hits", len(best_hits))
    
    best_hits = sorted(set(best_hits))

    with open(output_file, "w") as outfile:
        for hit in best_hits:
            outfile.write(hit + "\n")







class BiologicalSequence(ABC):
    def __init__(self, sequence:str):
        self._sequence = sequence

    def __len__(self):
        return len(self._sequence)
    
    def __getitem__(self, index):
        result = self._sequence[index]
        if isinstance(index, slice):
            return self.__class__(result)
        return result
    
    def __str__(self):
        return self._sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self._sequence}')"

    @abstractmethod
    def check_alphabet(self):
        ...

class NucleicAcidSequence(BiologicalSequence):
    _alphabet = None
    _complement_map = None  
    

    def check_alphabet(self):
       
        if self._alphabet is None:
            raise NotImplementedError("Alphabet not defined")
        return set(self._sequence).issubset(self._alphabet)

    def complement(self):
        if self._complement_map is None:
            raise NotImplementedError("Complement map not defined")
        return self.__class__(self._sequence.translate(self._complement_map))

    def reverse(self):
        return self.__class__(self._sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    _alphabet = set("ATGCatgc")
    _complement_map = str.maketrans("ATGCatgc", "TACGtacg")


    def transcribe(self):
        rna_seq = self._sequence.replace('T', 'U').replace('t', 'u')
        return RNASequence(rna_seq)

class RNASequence(NucleicAcidSequence):
    _alphabet = set("AUGCaugc")
    _complement_map = str.maketrans("AUGCaugc", "UACGuacg")


class AminoAcidSequence(BiologicalSequence):
    _alphabet = set("ACDEFGHIKLMNPQRSTVWY")


    def check_alphabet(self):
        return set(self._sequence).issubset(self._alphabet)

    def hydrophobic_ratio(self):
        hlen = len(self)
        hcnt = 0
        for s in self._sequence:
            if s in ("AVLIMFWP"):
                hcnt+=1
        return hcnt/hlen






def parse_blast():
    """
    Calls blast parsing function
    """
    input_file = "blast_results.txt"
    output_file = "best_hits.txt"

    parse_blast_results(input_file, output_file)
def extract_neighbor_genes(input_gbk: str, genes: str, output_fasta: str, n_before: int = 1, n_after: int = 1) -> str:


    """
    Extract protein sequences of neighboring genes from a GenBank file 
    and save them in FASTA format.
    """
    if isinstance(genes, str):
        genes = [genes]

    cds_list = []
    current_cds = {}
    inside_cds = False

    with open(input_gbk, "r") as f:
        for line in f:
            line = line.rstrip()

            if line.strip().startswith("CDS"):
                if current_cds:
                    cds_list.append(current_cds)
                current_cds = {}
                inside_cds = True

            elif inside_cds:
                if "/gene=" in line:
                    current_cds["gene"] = line.split("=")[1].strip('"')
                elif "/locus_tag=" in line:
                    current_cds["locus_tag"] = line.split("=")[1].strip('"')
                elif "/product=" in line:
                    current_cds["product"] = line.split("=")[1].strip('"')
                elif "/translation=" in line:
                    translation = line.split("=", 1)[1].strip().lstrip('"')
                    if translation.endswith('"'):
                        translation = translation.rstrip('"')
                        current_cds["translation"] = translation
                    else:
                        seq = [translation]
                        for next_line in f:
                            next_line = next_line.strip()
                            if next_line.endswith('"'):
                                seq.append(next_line.rstrip('"'))
                                break
                            else:
                                seq.append(next_line)
                        current_cds["translation"] = "".join(seq)

        if current_cds:
            cds_list.append(current_cds)


    extracted = []
    for i, cds in enumerate(cds_list):
        if cds.get("gene") in genes or cds.get("locus_tag") in genes:
            start = max(0, i - n_before)
            end = min(len(cds_list), i + n_after + 1)
            for j in range(start, end):
                if j == i:  
                    continue
                neighbor = cds_list[j]
                if "translation" in neighbor:
                    header = neighbor.get("locus_tag", neighbor.get("gene", "unknown"))
                    product = neighbor.get("product", "")
                    fasta_header = f">{header} {product}"
                    extracted.append((fasta_header, neighbor["translation"]))


    with open(output_fasta, "w") as out:
        for header, seq in extracted:
            out.write(header + "\n")
            out.write(seq + "\n")

def extract_neighbors():
    """
    Calls function to extract neighbors from .gbk
    """
    input_gbk = "ecoli.gbk"          
    genes_of_interest = ["blaTEM"] 
    output_fasta = "neighbors_only.fasta"

    extract_neighbor_genes(
        input_gbk=input_gbk,
        genes=genes_of_interest,
        n_before=2,
        n_after=2,
        output_fasta=output_fasta
    )

    print(f"Neighboring genes were extracted into {output_fasta}")








def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100),
                 length_bounds=(0, 2**32), quality_threshold=0):
    """
    Filters FASTQ reads using Biopython.
    """
    output_folder = "filtered"
    os.makedirs(output_folder, exist_ok=True)
    output_path = os.path.join(output_folder, output_fastq)

    good_reads = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        gc = gc_fraction(record.seq) * 100

        seq_len = len(record)

        qualities = record.letter_annotations["phred_quality"]
        mean_quality = sum(qualities) / len(qualities)


        if not (gc_bounds[0] <= gc <= gc_bounds[1]):
            continue
        if not (length_bounds[0] <= seq_len <= length_bounds[1]):
            continue
        if mean_quality < quality_threshold:
            continue

        good_reads.append(record)
    
    # LOGGING
    logger.info("Reads passing filters: %d", len(good_reads))

    SeqIO.write(good_reads, output_path, "fastq")



def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="biotools",
        description="Bioinformatics utilities: BLAST parsing, neighbor extraction, FASTQ filtering",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- blast ---
    p_blast = subparsers.add_parser("blast", help="Parse BLAST results and extract best hits")
    p_blast.add_argument("-i", "--input",  required=True, help="Input BLAST .txt file")
    p_blast.add_argument("-o", "--output", required=True, help="Output file for best hits")

    # --- neighbors ---
    p_nb = subparsers.add_parser("neighbors", help="Extract neighboring genes from a GenBank file")
    p_nb.add_argument("-i", "--input",  required=True, help="Input .gbk file")
    p_nb.add_argument("-g", "--genes",  required=True, nargs="+", help="Gene name(s) of interest")
    p_nb.add_argument("-o", "--output", required=True, help="Output .fasta file")
    p_nb.add_argument("--before", type=int, default=1, help="Neighbors to take before (default: 1)")
    p_nb.add_argument("--after",  type=int, default=1, help="Neighbors to take after (default: 1)")

    # --- filter ---
    p_fq = subparsers.add_parser("filter", help="Filter FASTQ reads by GC, length, and quality")
    p_fq.add_argument("-i", "--input",   required=True, help="Input .fastq file")
    p_fq.add_argument("-o", "--output",  required=True, help="Output .fastq filename (saved in filtered/)")
    p_fq.add_argument("--gc-min",        type=float, default=0,   help="Min GC%% (default: 0)")
    p_fq.add_argument("--gc-max",        type=float, default=100, help="Max GC%% (default: 100)")
    p_fq.add_argument("--len-min",       type=int,   default=0,   help="Min read length (default: 0)")
    p_fq.add_argument("--len-max",       type=int,   default=2**32, help="Max read length")
    p_fq.add_argument("--quality",       type=float, default=0,   help="Min mean Phred quality (default: 0)")

    return parser


if __name__ == "__main__":
    import logging

    logging.basicConfig(
        filename="biotools.log",
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    logger = logging.getLogger(__name__)

    parser = build_parser()
    args = parser.parse_args()

    try:
        if args.command == "blast":
            logger.info("Running BLAST parser: %s -> %s", args.input, args.output)
            parse_blast_results(args.input, args.output)
            logger.info("BLAST parsing complete. Output: %s", args.output)

        elif args.command == "neighbors":
            logger.info(
                "Extracting neighbors from %s for genes %s (before=%d, after=%d)",
                args.input, args.genes, args.before, args.after,
            )
            extract_neighbor_genes(
                input_gbk=args.input,
                genes=args.genes,
                output_fasta=args.output,
                n_before=args.before,
                n_after=args.after,
            )
            logger.info("Neighbor extraction complete. Output: %s", args.output)

        elif args.command == "filter":
            logger.info(
                "Filtering FASTQ %s (gc=%.1f-%.1f, len=%d-%d, quality>=%.1f)",
                args.input, args.gc_min, args.gc_max,
                args.len_min, args.len_max, args.quality,
            )
            filter_fastq(
                input_fastq=args.input,
                output_fastq=args.output,
                gc_bounds=(args.gc_min, args.gc_max),
                length_bounds=(args.len_min, args.len_max),
                quality_threshold=args.quality,
            )
            logger.info("FASTQ filtering complete. Output saved in filtered/%s", args.output)

    except FileNotFoundError as e:
        logger.error("File not found: %s", e)
        raise SystemExit(1)
    except Exception as e:
        logger.error("Unexpected error in command '%s': %s", args.command, e)
        raise SystemExit(1)