from src.helpers.fasta_tools import read_fasta

# Recreates identified mutations that improve performance identified through wet lab research

########################################################################################################################
# Single Amino acid mutations
########################################################################################################################

uvsx=read_fasta('../../data/uvsx/P04529_sequence.fasta')

# (E198N, E198R, E198K, K35G and K35G/E198R)
mutants=['E198N', 'E198R', 'E198K', 'K35G', 'K35G/E198R', 'D274A']

positions = {}

for mutant in mutants:
    parts = mutant.split('/')  # handles both single + double

    pos_list = []
    for part in parts:
        pos = int(part[1:-1]) - 1
        pos_list.append(pos)

    positions[mutant] = pos_list

for item in positions.items():
    new_base=item[0][-1]
    name=item[0]
    uvsx_list=list(uvsx)
    for pos in item[1]:
        uvsx_list[pos]=new_base
    new_seq=''.join(uvsx_list)
    safe_name = name.replace("/", "_")
    with open(f"../../data/uvsx/{safe_name}.fasta",'w') as f:
        f.write(f">{safe_name} Rational Design UvsX mutant\n{new_seq}")



########################################################################################################################
# Larger scale mutations
########################################################################################################################
# Getting sequence of loop 2 donor recombinase
# IN Future develop this to recreate published mutations involving larger scale modifications, in this case a full loop swap.
header_to_find = ">YP_003097304"
sequence_lines = []
found = False

with open('../data/ncbi/uvsx_sequences.fasta', 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if found:
                # We've reached the next header, stop reading
                break
            if line.startswith(header_to_find):
                found = True
        elif found:
            # Collect sequence lines
            sequence_lines.append(line)

# Combine all sequence lines into a single string
full_sequence = "".join(sequence_lines)
# print(f"{header_to_find}:\n{full_sequence}")
# idea for this was replicating the loop 2 swap from YP_003097304 into UvsX