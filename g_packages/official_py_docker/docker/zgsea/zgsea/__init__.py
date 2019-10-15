import pkg_resources
import json
import sys
import math
from collections import defaultdict
import pickle
import numpy as np
import pandas as pd
import scipy.stats
import colors
from tqdm import tqdm
from bisect import bisect_left
import statsmodels.stats.multitest


def progress(x, *args, **kwargs):
    print(x, *args, file=sys.stderr, **kwargs)


def info(x, *args, **kwargs):
    print(colors.blue(x), *args, file=sys.stderr, **kwargs)


def debug(x, *args, **kwargs):
    print(colors.green(x), *args, file=sys.stderr, **kwargs)
    pass


def parse_gmt(gmt_str):
    lines = gmt_str.split("\n")
    gsets = []
    for line in lines:
        if line == "":
            continue

        name, url, *gene_ids = line.split("\t")
        gsets.append({"name": name, "url": url, "gene_ids": gene_ids})
    return gsets


gset_group_id_to_filename = {
    "kegg": "c2.cp.kegg.v6.2.symbols.gmt",
    "go": "c5.all.v6.2.symbols.gmt",
}


def load_gsets(gset_group_id):
    data_fn = gset_group_id_to_filename[gset_group_id]
    gmt_str = pkg_resources.resource_string(__name__, f"data/{data_fn}")
    gmt_str = str(gmt_str, "utf-8")
    gsets = parse_gmt(gmt_str)

    return gsets


def squash(gene_scores, gamma=1, epsilon=1e-7):
    gene_scores = np.array(gene_scores)
    gene_scores = (gene_scores - np.min(gene_scores)) / (np.max(gene_scores) - np.min(gene_scores) + epsilon)
    gene_scores = np.pow(gene_scores, gamma)

    return gene_scores


def build_curve(gene_id_to_rank, genes_in_gset, sorted_genes):
    gene_ranks = sorted(set([gene_id_to_rank[id_] for id_ in genes_in_gset]))

    curve = []
    last = 0
    total_hit = 0
    total_miss = 0
    for r in gene_ranks:
        curve.append((r - last, abs(sorted_genes[r]["score"])))
        total_hit += abs(sorted_genes[r]["score"])
        total_miss += r - last
        last = r + 1
    if last <= len(sorted_genes):
        curve.append((len(sorted_genes) - last, 0))
        total_miss += len(sorted_genes) - last

    return curve, total_hit, total_miss


def calc_es_by_curve(curve, total_hit, total_miss):
    slope = total_hit / total_miss
    highest = 0
    lowest = 0
    dyy = 0
    x_high, y_high = 0, 0
    x_low, y_low = 0, 0
    for dx, dy in curve:
        x_high += dx
        y_high += dy

        high = y_high - x_high * slope
        if high > highest:
            highest = high

        x_low += dx
        y_low += dyy

        low = y_low - x_low * slope
        if low < lowest:
            lowest = low

        dyy = dy

    # TODO: add pos for leading-edge genes
    return highest / total_hit, -lowest / total_hit


def calc_es_by_definition(gene_id_to_rank, genes_in_gset, sorted_genes):
    p_hit = 0
    p_hits = [p_hit]
    p_miss = 0
    p_misses = [p_miss]

    for gene in sorted_genes:
        if gene["id_"] in genes_in_gset:
            p_hit += abs(gene["score"])
        else:
            p_miss += 1

        p_hits.append(p_hit)
        p_misses.append(p_miss)

    p_hits = np.array(p_hits, dtype=float)
    p_hits /= p_hit

    p_misses = np.array(p_misses, dtype=float)
    p_misses /= p_miss

    p_delta = p_hits - p_misses

    return np.max(np.abs(p_delta))


def to_human_homolog(gene_ids, species, type="symbol"):
    if species == "mouse":
        homolog_tb = pd.read_csv(pkg_resources.resource_stream(__name__, f"data/HOM_MouseHumanSequence.rpt"), sep='\t')
    else:
        raise ValueError()

    if type == "symbol":
        id_col = "Symbol"
    elif type == "entrez":
        id_col = "EntrezGene ID"
    else:
        raise ValueError()

    homo_id_col = "HomoloGene ID"

    gene_to_homo = homolog_tb.drop_duplicates(id_col).set_index(id_col)[homo_id_col]
    homo_to_human_gene = (
        homolog_tb[homolog_tb["NCBI Taxon ID"] == 10090].drop_duplicates(homo_id_col).set_index(homo_id_col)
    )

    homo_ids = gene_to_homo.reindex(gene_ids)
    human_gene_ids = homo_to_human_gene.reindex(homo_ids)

    return human_gene_ids


def simple_fisher(input_list, gset_group_id, background_list=None):
    gsets = load_gsets(gset_group_id)

    if background_list is None:
        progress('Generate background gene list from the gene sets ...')
        background_list = set.union(*[set(gset['gene_ids']) for gset in gsets])
        progress('Generate background gene list from the gene sets DONE')

    background_set = set(background_list)

    input_set = set(input_list).intersection(background_set)
    not_input_set = background_set.difference(input_set)

    result_rows = []
    for gset in tqdm(gsets):
        gset_genes_set = set(gset['gene_ids']).intersection(background_set)
        not_gset_genes_set = background_set.difference(gset_genes_set)

        tp = len(gset_genes_set & input_set)
        fp = len(gset_genes_set & not_input_set)
        fn = len(not_gset_genes_set & input_set)
        tn = len(not_gset_genes_set & not_input_set)
        odds_ratio, p_val = scipy.stats.fisher_exact([[tp, fp], [fn, tn]])
        result_rows.append(
            {
                'gene_set_name': gset['name'],
                'size': len(gset['gene_ids']),
                'p_val': p_val,
                'odds_ratio': odds_ratio,
                'n_overlaps': len(gset_genes_set & input_set),
                'overlapping_genes': ', '.join(gset_genes_set & input_set),
            }
        )

    result_tb = pd.DataFrame.from_records(result_rows)
    _, result_tb["fdr"] = statsmodels.stats.multitest.fdrcorrection(result_tb["p_val"])

    return result_tb


def gsea(gene_ids, gene_scores, gset_group_id, n_repeats=1000, min_n_genes_in_gset=10):
    sorted_genes = np.array(list(zip(gene_ids, gene_scores)), dtype=[("id_", object), ("score", float)])
    sorted_genes = np.sort(sorted_genes, order="score")
    gene_id_to_rank = dict(zip(sorted_genes["id_"], range(len(sorted_genes))))

    gsets = load_gsets(gset_group_id)

    # remove genes that are not in the input gene list
    gene_ids_set = set(gene_ids)
    gsets_filtered = []
    for gset in gsets:
        gset["gene_ids"] = sorted(gene_ids_set.intersection(set(gset["gene_ids"])))

        if len(gset["gene_ids"]) < min_n_genes_in_gset:
            # info(f"The gene set {gset['name']} contains too few genes in your gene list. Skip.")
            continue

        if gset["gene_ids"] == gene_ids_set:
            # info(f"The gene set {gset['name']} contains *all* genes in your gene list. Skip.")
            continue

        gsets_filtered.append(gset)
    gsets = gsets_filtered

    # group gsets by their sizes
    gsets_by_size = defaultdict(list)
    for gset in gsets:
        gsets_by_size[len(gset["gene_ids"])].append(gset)

    # calculate es and p_val by the gset size groups
    enriched_gsets = []
    for size, gsets_with_size in tqdm(gsets_by_size.items()):
        pos_list = []
        neg_list = []
        for repeat in range(n_repeats):
            genes_in_gset = np.random.choice(sorted_genes, size, replace=False)["id_"]
            curve, total_miss, total_hit = build_curve(gene_id_to_rank, genes_in_gset, sorted_genes)
            pos, neg = calc_es_by_curve(curve, total_miss, total_hit)
            pos_list.append(pos)
            neg_list.append(neg)
        mean_es_pos = sum(pos_list) / len(pos_list)
        mean_es_neg = sum(neg_list) / len(neg_list)
        sorted_nes_list = sorted(max(p / mean_es_pos, n / mean_es_neg) for p, n in zip(pos_list, neg_list))

        for gset in gsets_with_size:
            genes_in_gset = gset["gene_ids"]

            curve, total_miss, total_hit = build_curve(gene_id_to_rank, genes_in_gset, sorted_genes)
            pos, neg = calc_es_by_curve(curve, total_miss, total_hit)
            es = max(pos, neg)
            nes = max(pos / mean_es_pos, neg / mean_es_neg)

            p_val = 1 - bisect_left(sorted_nes_list, nes) / n_repeats

            enriched_gset = {
                "gset_name": gset["name"],
                "gset_size": len(gset["gene_ids"]),
                "genes": ", ".join(genes_in_gset),
                "es": es,
                "nes": nes,
                "mean_es_pos": mean_es_pos,
                "mean_es_neg": mean_es_neg,
                "p_val": p_val,
            }
            enriched_gsets.append(enriched_gset)

    if len(enriched_gsets) == 0:
        return None

    result_tb = pd.DataFrame.from_records(enriched_gsets).sort_values("p_val")
    _, result_tb["fdr"] = statsmodels.stats.multitest.fdrcorrection(result_tb["p_val"])

    return result_tb


def test():
    with open("/home/xzhu/Desktop/Marker score (Cluster 1 vs. rest) (2).json", "r") as f:
        data = (json.load(f)).items()
    # data = [
    #     # --- BACKGROUND
    #     ('ACVRL1', 0.5),
    #     ('ADCY3', 0.5),
    #     ('ALS2', 0.5),
    #     ('AQP3', 0.5),
    #     ('AQP4', 0.5),
    #     ('ATG12', 0.5),
    #     ('AVP', 0.5),
    #     ('AVPR2', 0.5),
    #     ('BCL2L11', 0.5),
    #     ('BIRC6', 0.5),
    #     ('BRAF', 0.5),
    #     ('CHP2', 0.5),
    #     ('CLTC', 0.5),
    #     ('COL2A1', 0.5),
    #     ('CREB3', 0.5),
    #     ('CREB3L1', 0.5),
    #     ('CREB3L2', 0.5),
    #     ('CREB3L3', 0.5),
    #     ('CREB3L4', 0.5),
    #     ('CSF2', 0.5),
    #     ('CYP1A1', 0.5),
    #     ('DCTN2', 0.5),
    #     ('DCTN5', 0.5),
    #     ('DSCAML1', 0.5),
    #     ('DVL3', 0.5),
    #     ('DYNLL1', 0.5),
    #     ('EGFL8', 0.5),
    #     ('FADD', 0.5),
    #     ('FCER1G', 0.5),
    #     ('FCGR3A', 0.5),
    #     ('FCGR3B', 0.5),
    #     ('GPC4', 0.5),
    #     ('GPC6', 0.5),
    #     ('GZMB', 0.5),
    #     ('HLA-G', 0.5),
    #     ('HOXB3', 0.5),
    #     ('HOXD1', 0.5),
    #     ('HOXD4', 0.5),
    #     ('HYAL1', 0.5),
    #     ('IFNA1', 0.5),
    #     ('IFNA17', 0.5),
    #     ('IKBKB', 0.5),
    #     ('IL12A', 0.5),
    #     ('IL12B', 0.5),
    #     ('IL8', 0.5),
    #     ('ISG15', 0.5),
    #     ('ITGB2', 0.5),
    #     ('JAG2', 0.5),
    #     ('KRT8', 0.5),
    #     ('MAP3K1', 0.5),
    #     ('MAP3K7', 0.5),
    #     ('MAPK10', 0.5),
    #     ('MAPK11', 0.5),
    #     ('MAPK13', 0.5),
    #     ('MAPK9', 0.5),
    #     ('MYCN', 0.5),
    #     ('NCR1', 0.5),
    #     ('NFAT5', 0.5),
    #     ('NRAS', 0.5),
    #     ('NSF', 0.5),
    #     ('PIK3R5', 0.5),
    #     ('PKD2', 0.5),
    #     ('PLK4', 0.5),
    #     ('POLR1A', 0.5),
    #     ('POLR1B', 0.5),
    #     ('POLR1D', 0.5),
    #     ('POLR1E', 0.5),
    #     ('POLR2A', 0.5),
    #     ('POLR2B', 0.5),
    #     ('POLR2C', 0.5),
    #     ('POLR2D', 0.5),
    #     ('POLR2J', 0.5),
    #     ('POLR2J2', 0.5),
    #     ('POLR2J3', 0.5),
    #     ('POLR2K', 0.5),
    #     ('POLR2L', 0.5),
    #     ('POLR3A', 0.5),
    #     ('POLR3B', 0.5),
    #     ('POLR3C', 0.5),
    #     ('POLR3D', 0.5),
    #     ('POLR3F', 0.5),
    #     ('POLR3G', 0.5),
    #     ('POLR3GL', 0.5),
    #     ('POLR3H', 0.5),
    #     ('POLR3K', 0.5),
    #     ('PRF1', 0.5),
    #     ('PRKX', 0.5),
    #     ('PSMA1', 0.5),
    #     ('PSMB9', 0.5),
    #     ('PSMC2', 0.5),
    #     ('PSMD1', 0.5),
    #     ('PSMD12', 0.5),
    #     ('PSMD13', 0.5),
    #     ('PSMF1', 0.5),
    #     ('PTPRR', 0.5),
    #     ('RAB11A', 0.5),
    #     ('RAB11B', 0.5),
    #     ('RAB5A', 0.5),
    #     ('RAB5B', 0.5),
    #     ('RAET1E', 0.5),
    #     ('RBP4', 0.5),
    #     ('RIPK1', 0.5),
    #     ('RSPO3', 0.5),
    #     ('SH2D1B', 0.5),
    #     ('SHC4', 0.5),
    #     ('SIKE1', 0.5),
    #     ('SOX6', 0.5),
    #     ('TBKBP1', 0.5),
    #     ('TBX1', 0.5),
    #     ('TIAM1', 0.5),
    #     ('UBC', 0.5),
    #     ('ULBP3', 0.5),
    #     ('VANGL1', 0.5),
    #     ('VAV1', 0.5),
    #     ('VAV2', 0.5),
    #     ('XRCC2', 0.5),
    #     # --- KEGG_STEROID_BIOSYNTHESIS
    #     ('SOAT1', 0.1),
    #     ('LSS', 0.1),
    #     ('SQLE', 0.1),
    #     ('EBP', 0.1),
    #     ('CYP51A1', 0.1),
    #     ('DHCR7', 0.1),
    #     ('CYP27B1', 0.1),
    #     ('DHCR24', 0.1),
    #     ('HSD17B7', 0.1),
    #     ('MSMO1', 0.1),
    #     ('FDFT1', 0.1),
    #     ('SC5DL', 0.1),
    #     ('LIPA', 0.1),
    #     # --- KEGG_TAURINE_AND_HYPOTAURINE_METABOLISM
    #     ('GGT7', 0.01),
    #     ('CSAD', 0.01),
    #     ('GGT5', 0.01),
    #     ('GGT1', 0.3),
    #     ('GAD1', 0.3),
    #     ('BAAT', 0.3),
    #     ('GAD2', 0.3),
    #     ('GGT6', 0.3),
    #     ('CDO1', 0.3),
    #     ('ADO', 0.3),
    # ]
    gene_ids, gene_scores = map(list, zip(*data))
    # gene_scores = np.random.rand(len(gene_ids))

    # print(pd.DataFrame({'id_': gene_ids, 'scores': gene_scores}).sort_values('scores'))

    print(gsea(gene_ids, gene_scores, "kegg"))
