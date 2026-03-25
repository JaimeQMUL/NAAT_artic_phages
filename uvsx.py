
from src.helpers.directory_creation import CreateDirectoryStructure
from src.extracting_references.extracting_uniprot_references import QuerySearchUniprot
from src.extracting_references.extracting_ncbi_references import QuerySearchNCBI
from src.extracting_references.extracting_pfam_domain_references import InterproSearchUniprot

protein_name='uvsx'


#Create directory structure to store all data associated for UvsX mining
CreateDirectoryStructure(protein_name)

# Use search terms to query Uniprot and save sequences found along with all avaliable metadata
QuerySearchUniprot("(recA OR uvsX OR recombinase) AND taxonomy_id:10239", protein_name)

# Use search terms to query NCBI databases, saving sequences along with all avaliable metadata
QuerySearchNCBI("(uvsx OR recA)[Protein] AND txid10239[Organism:exp]", protein_name)


### Add pfam identification step here

# Scanning Uniprot database for proteins annotated to contain the InterPro ids identified by Pfam
interpro_ids = ["IPR049428", "IPR049047"]
InterproSearchUniprot(interpro_ids, protein_name)

