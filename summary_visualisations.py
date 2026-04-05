
from src.visualisations.fasta_summary_plots import *
from src.helpers.data_wrangling import *

fig, ax =PlotSequenceLengths('uvsx/data/curated_database/cleaned_curated_database.fasta')
plt.show()

#Get taxonomy Lookups
ncbi_lookup=BuildTaxLookup('uvsx/data/curated_database/ncbi_query_search_metadata.csv')

#Get taxonomy lookups
uniprot_lookup=BuildTaxLookup('uvsx/data/curated_database/uniprot_query_search_metadata.json')

#Import current accessions
df = pd.read_csv('uvsx/data/curated_database/cleaned_metadata.csv')
print(f'Number of accessions: {len(df)}')

taxon_ids=GetTaxonIds(df, [ncbi_lookup, uniprot_lookup])
print(f'Number of taxon ids: {len(taxon_ids)}')

lineages=fetch_taxonomy_ranked_lineage_list(taxon_ids)

#Plot the Family breakdown
families=FindRanks(lineages, 'family')
family_counts=FindUniqueCounts(families)
fig, ax = PlotRankCounts(family_counts, rank='Family', top_n=10)
plt.show()


#Plot the Genus breakdown
genera=FindRanks(lineages, 'genus')
genus_counts=FindUniqueCounts(genera)
fig, ax = PlotRankCounts(genus_counts, rank='Genus', top_n=10)
plt.show()

#Plot the species breakdown
species=FindRanks(lineages, 'species')
species_counts=FindUniqueCounts(species)
fig, ax = PlotRankCounts(species_counts, rank='Species', top_n=10)
plt.show()

#Plot overall summary of
fig, ax= PlotRanksSummary(lineages)
plt.show()


