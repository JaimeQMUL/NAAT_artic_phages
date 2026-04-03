
import matplotlib.pyplot as plt



def PlotSequenceLenths(fasta_file):

    seq_lengths = []

    with open(fasta_file) as f:
        seq = ''
        header = None

        # Get sequence lengths
        for line in f:
            line = line.strip()

            if line.startswith('>'):
                # process previous sequence
                if seq:
                    seq_lengths.append(len(seq))


                header = line
                seq = ''
            else:
                seq += line

        # handle last sequence
        if seq:
            seq_lengths.append(len(seq))


    plt.hist(seq_lengths, bins=1000)
    plt.xlabel("Sequence length")
    plt.ylabel("Frequency")
    plt.title("Distribution of sequence lengths")
    plt.show()
