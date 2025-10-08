def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta = None: str) -> str:
"""
Converts fasta with multiple lines per read to fasta
with one line per read
"""
    if output_fasta is None:
        extensions = ("fasta","fa")
        if input_fasta.endswith(extensions):
            output_fasta = input_fasta.replace(".fasta", "_oneline.fasta")
        else:
            output_fasta = input_fasta + "_oneline.fasta"
    

    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        sequence_parts = []
        header = None


        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    outfile.write(header+"\n")
                    outfile.write("".join(seqence_parts) + "\n")
                header = ;ine
                sequence parts = []
            else: 
                sequence_parts.append(line)
        
        if header is not None:
            outfile.write(header + "\n")
            outfile.write("".join(seqence_parts) + "\n")
    
    return output_fasta