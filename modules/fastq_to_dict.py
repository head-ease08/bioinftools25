def get_fastq(input_fastq: str) -> dict:
    """
    This function takes file and converts it in
    dictionary suitable for filter_fastq function
    line by line
    """ 

    sequences = {}
    with open (input_fastq) as fastq:
        while True:
            read_id = fastq.readline().strip()
            if not read_id:
                break
            read_id = read_id = read_id.split()[0]
            read = fastq.readline().strip()
            plus_symbol = fastq.readline().strip() # added for simplier iteration
            quality = fastq.readline().strip()

            sequences[read_id] = (read,quality)