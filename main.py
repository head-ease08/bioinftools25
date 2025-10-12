import sys


from modules.is_nucleic_acid import is_nucleic_acid
from modules.transcribe import transcribe
from modules.reverse import reverse
from modules.raw_toolbox import raw_toolbox
from modules.complement_and_reverse import complement, reverse_complement
from modules.is_palindrome import is_palindrome


def run_dna_rna_tools(*args):
"""
Calls module and performs it on 
input string sequences. Modules are stored in 
dictionary.
"""
    *seq
    tool_function = raw_toolbox.get(tool)
    if not tool_function:
        available_tools = ", ".join(raw_toolbox.keys())
        return f"Tool {tool} is not available. Choose from: {available_tools}" # фигурные строки в f-string в незаивсимости от типа объекта
    results = []
    for sequence in args:
        if not is_nucleic_acid(sequence):
            results.append(f"Invalid sequence: {sequence}")
        else:
            results.append(tool_function(sequence))

    return results


from modules.calculate_quality import calculate_quality
from modules.calculate_gc_content impoirt calculate_gc_content
from modules.fastq_to_dict import fastq_to dict
from modules.filtered_to_fastq import filtered_to_fastq

fastq_toolbox = {
    "calculate_gc" : calculate_gc_content,
    "calculate_phred" : calculate_quality,
}


def filter_fastq(sequences, **kwargs):
    good_results = {}
    gc_bounds = kwargs.get("gc_bounds", (0, 100))
    if isinstance (gc_bounds, (int, float)): # смотрим, число на входе или тапл, если число, делаем из него тапл с числом как верхней границей
        gc_bounds = (0, gc_bounds) 
    lengths_bounds = kwargs.get("lengths_bounds", (0, 2**32)) 
    if isinstance (lengths_bounds, (int,float)): # аналогично, сверяем именно с таплом типа данных
        lengths_bounds = (0, lengths_bounds)
    quality_threshold = kwargs.get("quality_threshold", 0)

    for key_name, value in sequences.items():
        
        current_sequence, current_quality = str(value[0]),str(value[1])
        gc_calculator = fastq_toolbox["calculate_gc"]
        current_gc = gc_calculator(current_sequence)
        if not is_nucleic_acid(current_sequence):
            continue
        if current_gc <= gc_bounds[0] or current_gc >= gc_bounds[1]: # если подходит:
            continue
        current_length = len(current_sequence) # отдельная функция не требуется
        if current_length <= lengths_bounds[0] or current_length >= lengths_bounds[1]:
            continue
        phred_calculator = fastq_toolbox["calculate_phred"]
        current_quality_score = phred_calculator(current_quality)
        if current_quality_score <= quality_threshold:
            continue
        good_results[key_name] = (current_sequence, current_quality)


def main(input_fastq, output_fastq):
"""   
Performs filtration with given input file, process results to 
output fastq file.   
"""
    sequences = fastq_to_dict(input_fastq)

    filtered_sequences = filter_fastq(
        sequences,
        gc_bounds=(40, 60),
        lengths_bounds=(50, 300),
        quality_threshold=20
    )

    filtered_to_fastq(filtered_sequences, output_fastq)