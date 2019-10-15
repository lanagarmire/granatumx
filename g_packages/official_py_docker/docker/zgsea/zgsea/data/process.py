import pandas as pd
print()

def mouse_to_human_homolog(homolog_tb, gene_ids, gene_scores, type='symbol'):
    possible_ids = homolog_tb.loc[homolog_tb['Symbol'] == mouse_gene_id, 'HomoloGene ID']
    if len(possible_ids) > 0:
        homo_id = possible_ids[0]
    else:
        return None
    if type == symbol

df = pd.read_csv('./HOM_MouseHumanSequence.rpt', sep='\t')
for k, gdf in df.groupby('HomoloGene ID'):
    if sum(gdf['NCBI Taxon ID'] == 9606) > 1:
        print(gdf)
        break
