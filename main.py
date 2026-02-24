import sys
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os




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

    best_hits = sorted(set(best_hits))

    with open(output_file, "w") as outfile:
        for hit in best_hits:
            outfile.write(hit + "\n")





from abc import ABC, abstractmethod

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


    SeqIO.write(good_reads, output_path, "fastq")