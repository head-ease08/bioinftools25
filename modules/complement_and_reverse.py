complement_map = str.maketrans("ATGCatgcUu", "TACGtacgAa")


def complement(sequence: str) -> str:
    """
    Reads string nucleic acid sequence and returns complement nucleic acid string. 
    Complement string is built with translation.
    """
    return sequence.translate(complement_map)




def reverse_complement(sequence: str) -> str:
    """
    Reads string nucleic acid sequence and returns reversed complement nucleic acid string. 
    Complement string is built with translation.
    """
    return reverse(complement(sequence))

