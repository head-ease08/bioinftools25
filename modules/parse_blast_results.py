def parse_blast_results(input_file: str, output_file: str) -> str:
    """
    Reads BLAST txt results and extracts the best hit (first in alignments block)
    for each query. Saves unique protein descriptions in one-column file.
    """
    best_hits = []

    with open(input_file, "r") as infile:
        lines = infile.readlines()

    for i, line in enumerate(lines):
        if line.strip().startswith("Sequences producing significant alignments:"):
            if i + 1 < len(lines):
                first_hit = lines[i + 1].strip()
                if first_hit:
                    parts = first_hit.split(maxsplit=1)
                    if len(parts) > 1:
                        description = parts[1]
                    else:
                        description = parts[0]
                    best_hits.append(description)

    best_hits = sorted(set(best_hits))

    with open(output_file, "w") as outfile:
        for hit in best_hits:
            outfile.write(hit + "\n")