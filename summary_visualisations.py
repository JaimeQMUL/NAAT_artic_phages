import matplotlib.pyplot as plt

from src.visualisations.fasta_summary_plots import *
from src.helpers.data_wrangling import *



#Get ncbi taxonomy Lookups
ncbi_lookup=BuildTaxLookup('uvsx/data/curated_database/ncbi_query_search_metadata.csv')

#Get uniprot taxonomy lookups
uniprot_lookup=BuildTaxLookup('uvsx/data/curated_database/uniprot_query_search_metadata.json')

########################################################################################################################
# Plots for cleaned database
########################################################################################################################
fig, ax = PlotSequenceLengths('uvsx/data/curated_database/cleaned_curated_database.fasta')
plt.savefig('uvsx/plots/cleaned_sequence_lengths.png', dpi=300, bbox_inches='tight')


#Import current accessions
df = pd.read_csv('uvsx/data/curated_database/cleaned_metadata.csv')
print(f'Number of accessions: {len(df)}')

taxon_ids=GetTaxonIds(df, [ncbi_lookup, uniprot_lookup])
print(f'Number of taxon ids: {len(taxon_ids)}')

lineages=fetch_taxonomy_ranked_lineage_list(taxon_ids)

#Plot the Family breakdown
families=FindRanks(lineages, 'family')
family_counts=FindUniqueCounts(families)
fam_fig, ax = PlotRankCounts(family_counts, rank='Family', top_n=10)
fam_fig.savefig("uvsx/plots/cleaned_family_counts.png", dpi=300, bbox_inches='tight')


#Plot the Genus breakdown
genera=FindRanks(lineages, 'genus')
genus_counts=FindUniqueCounts(genera)
gen_fig, ax = PlotRankCounts(genus_counts, rank='Genus', top_n=10)
gen_fig.savefig("uvsx/plots/cleaned_genus_counts.png", dpi=300, bbox_inches='tight')


#Plot the species breakdown
species=FindRanks(lineages, 'species')
species_counts=FindUniqueCounts(species)
spec_fig, ax = PlotRankCounts(species_counts, rank='Species', top_n=10)
spec_fig.savefig("uvsx/plots/cleaned_species_counts.png", dpi=300, bbox_inches='tight')

#Plot overall summary of
sum_fig, ax= PlotRanksSummary(lineages)
sum_fig.savefig("uvsx/plots/cleaned_ranks_summary.png", dpi=300, bbox_inches='tight')


########################################################################################################################
# Plots for filtered database
########################################################################################################################
fig, ax = PlotSequenceLengths('uvsx/data/curated_database/filtered_curated_database.fasta')
plt.savefig('uvsx/plots/filtered_sequence_lengths.png', dpi=300, bbox_inches='tight')


#Import current accessions
df = pd.read_csv('uvsx/data/curated_database/filtered_metadata.csv')
print(f'Number of accessions: {len(df)}')

taxon_ids=GetTaxonIds(df, [ncbi_lookup, uniprot_lookup])
print(f'Number of taxon ids: {len(taxon_ids)}')

lineages=fetch_taxonomy_ranked_lineage_list(taxon_ids)

#Plot the Family breakdown
families=FindRanks(lineages, 'family')
family_counts=FindUniqueCounts(families)
fam_fig, ax = PlotRankCounts(family_counts, rank='Family', top_n=10)
fam_fig.savefig("uvsx/plots/filtered_family_counts.png", dpi=300, bbox_inches='tight')


#Plot the Genus breakdown
genera=FindRanks(lineages, 'genus')
genus_counts=FindUniqueCounts(genera)
gen_fig, ax = PlotRankCounts(genus_counts, rank='Genus', top_n=10)
gen_fig.savefig("uvsx/plots/filtered_genus_counts.png", dpi=300, bbox_inches='tight')


#Plot the species breakdown
species=FindRanks(lineages, 'species')
species_counts=FindUniqueCounts(species)
spec_fig, ax = PlotRankCounts(species_counts, rank='Species', top_n=10)
spec_fig.savefig("uvsx/plots/filtered_species_counts.png", dpi=300, bbox_inches='tight')

#Plot overall summary of
sum_fig, ax= PlotRanksSummary(lineages)
sum_fig.savefig("uvsx/plots/filtered_ranks_summary.png", dpi=300, bbox_inches='tight')




########################################################################################################################
# Plots for top hits database
########################################################################################################################
fig, ax = PlotSequenceLengths('uvsx/data/curated_database/top_hits.fasta')
plt.savefig('uvsx/plots/top_hits_sequence_lengths.png', dpi=300, bbox_inches='tight')


#Import current accessions
df = pd.read_csv('uvsx/data/curated_database/top_hits.csv')
print(f'Number of accessions: {len(df)}')

taxon_ids=GetTaxonIds(df, [ncbi_lookup, uniprot_lookup])
print(f'Number of taxon ids: {len(taxon_ids)}')

lineages=fetch_taxonomy_ranked_lineage_list(taxon_ids)

#Plot the Family breakdown
families=FindRanks(lineages, 'family')
family_counts=FindUniqueCounts(families)
fam_fig, ax = PlotRankCounts(family_counts, rank='Family', top_n=10)
fam_fig.savefig("uvsx/plots/top_hits_family_counts.png", dpi=300, bbox_inches='tight')


#Plot the Genus breakdown
genera=FindRanks(lineages, 'genus')
genus_counts=FindUniqueCounts(genera)
gen_fig, ax = PlotRankCounts(genus_counts, rank='Genus', top_n=10)
gen_fig.savefig("uvsx/plots/top_hits_genus_counts.png", dpi=300, bbox_inches='tight')


#Plot the species breakdown
species=FindRanks(lineages, 'species')
species_counts=FindUniqueCounts(species)
spec_fig, ax = PlotRankCounts(species_counts, rank='Species', top_n=10)
spec_fig.savefig("uvsx/plots/top_hits_species_counts.png", dpi=300, bbox_inches='tight')

#Plot overall summary of
sum_fig, ax= PlotRanksSummary(lineages)
sum_fig.savefig("uvsx/plots/top_hits_ranks_summary.png", dpi=300, bbox_inches='tight')

