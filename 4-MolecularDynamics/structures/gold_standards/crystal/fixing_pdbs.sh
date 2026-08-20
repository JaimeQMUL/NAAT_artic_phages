
# Fixing crystal structures pdbs
pdbfixer 3IO5.pdb --output 3IO5.fix.pdb
pdbfixer 7Z3M.pdb --output 7Z3M.fix.pdb
pdbfixer 9GBG.pdb --output 9GBG.fix.pdb

python3 <<'EOF'
from modeller import *
from modeller.scripts import complete_pdb

env = Environ()

env.libs.topology.read(file="/opt/anaconda3/lib/modeller-10.7/modlib/top_heav.lib")
env.libs.parameters.read(file="/opt/anaconda3/lib/modeller-10.7/modlib/par.lib")

mdl = complete_pdb(env, "7Z3M.pdb")

mdl.write(file="7Z3M.complete.pdb")

print("Created 7Z3M.complete.pdb")
EOF


python3 - <<'EOF'
inp = "7Z3M.complete.pdb"
out = "7Z3M.fixedformat.pdb"

with open(inp) as f, open(out, "w") as o:
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            resfield = line[17:22].strip()   # e.g. "HISA" or "HISA "
            resname = resfield[:3]
            chain = resfield[3] if len(resfield) > 3 else "A"
            resid = line[22:26].strip()

            newline = (
                line[:17] +
                f"{resname:>3}" +   # cols 18-20
                " " +               # col 21 (blank)
                f"{chain:1}" +      # col 22 chainID
                f"{int(resid):4d}" +# cols 23-26 resSeq
                line[26:]
            )
            o.write(newline)
        else:
            o.write(line)

print("done")
EOF


# Removing second dimer from 3IO5
awk 'BEGIN{keep=1} /^ATOM|^HETATM|^ANISOU/ {keep=(substr($0,22,1)=="A")} keep' 3IO5.fix.pdb > 3IO5_A.pdb
echo "END" >> 3IO5_A.pdb



