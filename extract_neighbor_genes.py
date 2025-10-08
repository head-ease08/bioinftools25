def extract_neighbor_genes(input_gbk, genes, output_fasta, n_before=1, n_after=1):
    if isinstance(genes, str):
        genes = [genes]

    cds_list = []
    current_cds = {}
    inside_cds = False

    with open(input_gbk, "r") as f:
        for line in f:
            line = line.rstrip()

            if line.strip().startswith("CDS"):
                if current_cds:
                    cds_list.append(current_cds)
                current_cds = {}
                inside_cds = True

            elif inside_cds:
                if "/gene=" in line:
                    current_cds["gene"] = line.split("=")[1].strip('"')
                elif "/locus_tag=" in line:
                    current_cds["locus_tag"] = line.split("=")[1].strip('"')
                elif "/product=" in line:
                    current_cds["product"] = line.split("=")[1].strip('"')
                elif "/translation=" in line:
                    translation = line.split("=", 1)[1].strip().lstrip('"')
                    if translation.endswith('"'):
                        translation = translation.rstrip('"')
                        current_cds["translation"] = translation
                    else:
                        seq = [translation]
                        for next_line in f:
                            next_line = next_line.strip()
                            if next_line.endswith('"'):
                                seq.append(next_line.rstrip('"'))
                                break
                            else:
                                seq.append(next_line)
                        current_cds["translation"] = "".join(seq)

        if current_cds:
            cds_list.append(current_cds)


    extracted = []
    for i, cds in enumerate(cds_list):
        if cds.get("gene") in genes or cds.get("locus_tag") in genes:
            start = max(0, i - n_before)
            end = min(len(cds_list), i + n_after + 1)
            for j in range(start, end):
                if j == i:  
                    continue
                neighbor = cds_list[j]
                if "translation" in neighbor:
                    header = neighbor.get("locus_tag", neighbor.get("gene", "unknown"))
                    product = neighbor.get("product", "")
                    fasta_header = f">{header} {product}"
                    extracted.append((fasta_header, neighbor["translation"]))


    with open(output_fasta, "w") as out:
        for header, seq in extracted:
            out.write(header + "\n")
            out.write(seq + "\n")