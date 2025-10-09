import os


def filtered_to_fastq(sequences: str, output_fastq:str) -> str:
    """
    Writes filtered reads to FASTQ file.
    """
    output_folder = "filtered"
    os.makedirs(output_folder, exist_ok = True)
    output_path = os.path.join(output_folder, output_fastq)

    if os.path.exists(output_path):
        raise FileExistsError(f"File {output_path} exists")
    
    with open(output_path, "w") as f:
        for read_id, (read,quality) in sequences.items():
            f.write(f"{read_id}\n")
            f.write(f"{read}\n")
            f.write("+\n")
            f.write(f"{quality}\n")