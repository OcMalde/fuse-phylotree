#!/usr/bin/env python3

import argparse
import hdbscan
import pandas as pd
import numpy as np
from kmodes.kmodes import KModes
import matplotlib.pyplot as plt
import scipy as scipy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist
from sklearn.metrics.cluster import adjusted_rand_score

#==============================================================================
# Load datas
#==============================================================================

def load_fchange_mchange_expand_csv(filename):
    """
    """
    dict_gene_modulesGained = {}
    with open(filename, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n","").replace(" ","").split(",")
            if "gene" in splited_line: continue
            modules_gained = [mod for mod in splited_line[4].split("|") if mod != ""]
            if len(modules_gained) > 0:
                dict_gene_modulesGained[str(splited_line[0] + " " + splited_line[1])] = modules_gained
    return dict_gene_modulesGained

def build_panda_dataframe(dict_gene_modulesList):
    """
    """
    gene_list = []
    all_modules_list = []
    for gene, modules_list in dict_gene_modulesList.items():
        for module in modules_list:
            if module not in all_modules_list:
                all_modules_list.append(module)
        gene_list.append(gene)
    dict_mod_array = {}
    for module in all_modules_list:
        mod_list = []
        for gene in gene_list:
            if module in dict_gene_modulesList[gene]:
                mod_list.append("1")
            else:
                mod_list.append("0")
        dict_mod_array[module] = mod_list
    dict_mod_array["gene"] = gene_list
    data = pd.DataFrame(dict_mod_array)
    data = data.set_index("gene")
    return data

#==============================================================================
# Find optimal K value
#==============================================================================

def elbow_curve(data):
    """
    https://www.analyticsvidhya.com/blog/2021/06/kmodes-clustering-algorithm-for-categorical-data/
    """
    # Elbow curve to find optimal K
    cost = []
    maxi = len(data.index)
    maxi = len([gene_leaf.split(" ")[0] for gene_leaf in data.index])
    print(maxi)
    K = range(1, maxi)
    for num_clusters in list(K):
        kmode = KModes(n_clusters=num_clusters, init = "random", n_init = 50, verbose=1)
        kmode.fit_predict(data)
        cost.append(kmode.cost_)
    plt.plot(K, cost, 'bx-')
    plt.xlabel('No. of clusters')
    plt.ylabel('Cost')
    plt.title('Elbow Method For Optimal k')
    plt.savefig("elbow.pdf", dpi=300, format="pdf")
    plt.show()

#==============================================================================
# Build clusters
#==============================================================================

def kmodes_k_clusters(K, data):
    """
    """
    kmode = KModes(n_clusters=K, init = "random", n_init = 5, verbose=1)
    clusters = kmode.fit_predict(data)
    data.insert(0, "Cluster", clusters, True)
    ari(data)
    return data

def hdbscan_clusters(data):
    """
    https://hdbscan.readthedocs.io/en/latest/basic_hdbscan.html
    """
    clusterer = hdbscan.HDBSCAN(metric='hamming')
    clusterer.fit(data)
    data.insert(0, "Cluster", clusterer.labels_, True)
    ari(data)
    return data

#==============================================================================
# Adjusted Rand Index (compare clusters)
#==============================================================================

def ari(data):
    """
    https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html
    https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index
    """
    gene_clust = [gene_leaf.split(" ")[0] for gene_leaf in data.index] 
    clust = data["Cluster"]
    ari_score = adjusted_rand_score(gene_clust,clust)
    print(f"Adjusted Rand Index score (ARI) : {ari_score} (0-1)")

#==============================================================================
# SLINK
#==============================================================================

def slink(data):
    """
    https://stackoverflow.com/questions/60016770/separating-clusters-after-using-slink-in-python-r
    """
    N = 150
    dist = pdist(data, metric="hamming")
    Z = linkage(data, "single")
    plt.figure(figsize=(1000,1000))
    plt.figure()
    plt.xticks(range(N)) # add loads of ticks
    plt.grid()
    plt.gca().margins(x=0)
    plt.gcf().canvas.draw()
    tl = plt.gca().get_xticklabels()
    maxsize = max([t.get_window_extent().width for t in tl])
    m = 0.8 # inch margin
    s = maxsize/plt.gcf().dpi*N+2*m
    margin = m/plt.gcf().get_size_inches()[0]
    plt.gcf().subplots_adjust(left=margin, right=1.-margin)
    plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])
    dn = scipy.cluster.hierarchy.dendrogram(Z, labels=[g.split("_")[0] for g in data.index], orientation="top", leaf_font_size=7, leaf_rotation=90)
    plt.tight_layout()
    plt.savefig("slink.pdf", dpi=300, format="pdf")

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("modules_change_expand",
                        help = "csv file containing modules compositions changes expanded to the leaves modules composition for each internal node (ex : functionChange_moduleChange_expand_seadogMD_712.csv)",
                        type=str)
    parser.add_argument("--elbow_curve",
                        help = "Draw the elbow curve to search for an optimal K",
                        action="store_true")
    parser.add_argument("-K",
                        help = "K (i.e.; number of Kmodes clusters to search)",
                        type=int)
    parser.add_argument("--HDBSCAN",
                        help = "HDBSCAN clustering",
                        action="store_true")
    parser.add_argument("--SLINK",
                        help = "Single Link Hierarchical clustering",
                        action="store_true")
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

if __name__ == "__main__":
    args = parser()
    dict_gene_modulesGained = load_fchange_mchange_expand_csv(args.modules_change_expand)
    data = build_panda_dataframe(dict_gene_modulesGained)
    if args.elbow_curve:
        elbow_curve(data)
    if args.K:
        K = args.K
        data = kmodes_k_clusters(K, data)
        data.to_csv(f"kmodes_K{K}.csv")
    if args.HDBSCAN:
        data = hdbscan_clusters(data)
        data.to_csv(f"hdbscan_clust.csv")
    if args.SLINK:
        slink(data)


