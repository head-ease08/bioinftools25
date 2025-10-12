def transcribe(sequence:str) -> str:
    """
    Reads string DNA sequence and returns transcribed RNA string.
    Reads string RNA sequence and returns transcribed RNA string.
    """
    return sequence.replace("t", "u").replace("T", "U")

