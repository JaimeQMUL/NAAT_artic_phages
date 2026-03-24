from pathlib import Path


def read_fasta(filename):

    path_to_data = Path(filename)
    seq=''
    with open(path_to_data, 'r') as f:
        for l in f.readlines():
            if l.startswith('>') == False:
                seq += l.rstrip().upper()
    return seq