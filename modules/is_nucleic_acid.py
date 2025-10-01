def is_nucleic_acid(sequence: str) -> bool:
"""
Reads input string and checks if it consists only from 
allowed nucleotides (stored in 2 sets for DNA and RNA
respectively).
Returns boolean value.
"""
    dna_nucleotides = set("AGCTagct")
    rna_nucleotides = set("ACTactuU")
    input_set = set(sequence)
    return input_set <= (dna_nucleotides) or input_set <= (rna_nucleotides)

