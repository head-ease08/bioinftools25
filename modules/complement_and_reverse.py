complement_map = str.maketrans("ATGCatgcUu", "TACGtacgAa")


def complement(sequence: str) -> str:
    return sequence.translate(complement_map)

"""
Reads string nucleic acid sequence and returns complement nucleic acid string. 
Complement string is built with translation.
"""


def reverse_complement(sequence: str) -> str:
    return reverse(complement(sequence))

"""
Reads string nucleic acid sequence and returns reversed complement nucleic acid string. 
Complement string is built with translation.
"""