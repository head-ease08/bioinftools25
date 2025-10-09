def is_palindrome(sequence: str) -> bool:
    """
    Checks if string nucleic acid sequence is palindrome or not.
    Returns boolean value.
    """
    return bool(sequence.upper() == reverse_complement(sequence).upper())

