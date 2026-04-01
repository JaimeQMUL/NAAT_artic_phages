from pathlib import Path


def read_fasta(filename):

    path_to_data = Path(filename)
    seq=''
    with open(path_to_data, 'r') as f:
        for l in f.readlines():
            if l.startswith('>') == False:
                seq += l.rstrip().upper()
    return seq


# Create a function that takes an accession and database file, and returns the sequence
def ExtractSequence(accession, filename, header=False):
    seq = ""
    capture = False
    matched_header = None

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                current_header = line
                first_half = len(line) // 2
                check_header=current_header[:first_half]


                # Extract accession from UniProt/NCBI header formats

                if accession in check_header:
                    print(f'Accession: {accession} Header: {current_header}')
                    capture = True
                    matched_header = current_header
                    continue
                else:
                    if capture:
                        break
                    capture = False

            elif capture:
                seq += line

    if header:
        return matched_header, seq

    return seq


