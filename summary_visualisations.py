
from src.visualisations.fasta_summary_plots import *
from src.helpers.data_wrangling import *

fig, ax = PlotSequenceLengths('uvsx/data/curated_database/cleaned_curated_database.fasta')
plt.show()

#Get ncbi taxonomy Lookups
ncbi_lookup=BuildTaxLookup('uvsx/data/curated_database/ncbi_query_search_metadata.csv')

#Get uniprot taxonomy lookups
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
fam_fig, ax = PlotRankCounts(family_counts, rank='Family', top_n=10)
fam_fig.savefig("uvsx/plots/family_counts.png", dpi=300, bbox_inches='tight')



#Plot the Genus breakdown
genera=FindRanks(lineages, 'genus')
genus_counts=FindUniqueCounts(genera)
gen_fig, ax = PlotRankCounts(genus_counts, rank='Genus', top_n=10)
gen_fig.savefig("uvsx/plots/genus_counts.png", dpi=300, bbox_inches='tight')


#Plot the species breakdown
species=FindRanks(lineages, 'species')
species_counts=FindUniqueCounts(species)
spec_fig, ax = PlotRankCounts(species_counts, rank='Species', top_n=10)
spec_fig.savefig("uvsx/plots/species_counts.png", dpi=300, bbox_inches='tight')

#Plot overall summary of
sum_fig, ax= PlotRanksSummary(lineages)
sum_fig.savefig("uvsx/plots/ranks_summary.png", dpi=300, bbox_inches='tight')


