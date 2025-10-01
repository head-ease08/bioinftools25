def calculate_gc_content(sequence: str) -> float:
    """
    Calculates content of guanine and cytosine in nucleic acid string.
    If string is empty returns zero.
    Takes string on input, returns float
    """
    if not sequence: 
        return 0.0
    
    gc_count = 0
    for nucleotide in sequence.upper():
        if nucleotide in ('G', 'C'):
            gc_count += 1
            
    return (gc_count / len(sequence)) * 100