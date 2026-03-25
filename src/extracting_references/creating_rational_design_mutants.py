from src.helpers.fasta_tools import read_fasta

# Recreates identified mutations that improve performance identified through wet lab research

uvsx=read_fasta('../../data/uvsx/P04529_sequence.fasta')

# (E198N, E198R, E198K, K35G and K35G/E198R)
mutants=['E198N', 'E198R', 'E198K', 'K35G', 'K35G/E198R']

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
        f.write(f">{safe_name} Rational Design UvsX mutant https://pubs.acs.org/doi/10.1021/acs.biochem.5c00098\n{new_seq}")


