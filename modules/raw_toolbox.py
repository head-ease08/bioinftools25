from modules.is_nucleic_acid import is_nucleic_acid
from modules.transcribe import transcribe
from modules.reverse import reverse
from modules.complement_and_reverse import complement, reverse_complement
from modules.is_palindrome import is_palindrome

raw_toolbox = {
    """
    Dictionary used to call functions in main.py
    """
    "is_na": is_nucleic_acid,
    "transcribe": transcribe,
    "reverse": reverse,
    "reverse_complement": reverse_complement,
    "is_palindrome": is_palindrome
}

