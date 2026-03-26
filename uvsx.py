
from src.helpers.directory_creation import CreateDirectoryStructure
from src.extracting_references.creating_curated_database import *
from src.biotools.fasta_tools import read_fasta
from src.helpers.uniprot_api import stage_uniprot_protein


protein_name='uvsx'
accession='P04529'

# Get reference protein sequence
stage_uniprot_protein(accession, f'{protein_name}/data/references')

#Create directory structure to store all data associated for UvsX mining
CreateDirectoryStructure(protein_name)

# Use search terms to query Uniprot and save sequences found along with all avaliable metadata
QuerySearchUniprot("(recA OR uvsX OR recombinase) AND taxonomy_id:10239", protein_name)

sleep(2)

# Use search terms to query NCBI databases, saving sequences along with all avaliable metadata
QuerySearchNCBI("(uvsx OR recA)[Protein] AND txid10239[Organism:exp]", protein_name)


### Add pfam identification step here. use model to pull out domain IDs. Or get them from uniprot/ncbi.

# Scanning Uniprot database for proteins annotated to contain the InterPro ids identified by Pfam
interpro_ids = ["IPR049428", "IPR049047"]
InterproSearchUniprot(interpro_ids, protein_name)
# I find 1048 instead of 1302

# From identified mutants in the literature, creating sequences with the mutations identified.
mutants=['E198N', 'E198R', 'E198K', 'K35G', 'K35G/E198R', 'D274A'] #mutants identified from literature
reference=read_fasta(f'{protein_name}/data/references/P04529_sequence.fasta') #getting protein sequence of UvsX reference
CreateMutants(mutants, reference, protein_name) #Creating reference sequences and writing them to a fasta file


# Extracting mutants identified from literature and stored in the RCSB database
# Sequences identified in paper from ice samples. 7Z3M is thermophilic and 9GBG is psychrophilic  https://academic.oup.com/nar/article/54/4/gkag069/8462617?login=true
pdb_ids=['7Z3M','9GBG']
IDSearchRCSB(pdb_ids, protein_name)