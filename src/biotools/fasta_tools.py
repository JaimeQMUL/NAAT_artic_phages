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
                    # print(f'Accession: {accession} Header: {current_header}') # debug to check accession matches header
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

def CreateDomainFasta(results_file, output_file, reference_file):

    id=output_file.split('/')[-1].split('_')[0]

    coords = []
    with open(results_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()

            seq_id = fields[3]
            evalue = float(fields[12])
            start = int(fields[17])
            end = int(fields[18])

            if evalue < 1e-5:
                coords.append((seq_id, start, end))

        print(coords)


    with open(output_file, 'w') as f:
        for c in coords:
            start=c[1]-1
            end=c[2]
            seq=read_fasta(reference_file).strip()
            domain=seq[start:end]

            f.write(f'>{c[0]} {id} Domain \n{domain}\n\n')

