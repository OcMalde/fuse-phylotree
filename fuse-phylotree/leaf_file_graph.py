#!/bin/python3

# Make a network x graph of protein-protein interactions from a psicquic

import argparse
from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt


def parse_ppi_csv(ppi_file) -> dict:
    """
    Read and parse a leaf csv file
    """
    dict_prot_phenotypeList = {}
    with open(ppi_file, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n","").split(",")
            prot = splited_line[0] 
            phenotype_list = splited_line[1].split("|")
            dict_prot_phenotypeList[prot] = phenotype_list
    return dict_prot_phenotypeList

def load_equivalence(equivalence_csv) -> dict:
    """
    Load equivalence csv file, 
    return a dict 
    { name : refseq }
    """
    dict_name_refseq = {}
    with open(equivalence_csv, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n","").split(",")
            name = splited_line[0]
            refseq = splited_line[1]
            dict_name_refseq[name] = refseq
    return dict_name_refseq

def make_interaction_list(dict_prot_phenotypeList, dict_name_refseq) -> list:
    raw_interaction_list = []
    interaction_list = []
    all_phenotype_list = []
    for prot, phenotype_list in dict_prot_phenotypeList.items():
        for name, refseq in dict_name_refseq.items():
            if prot == refseq: prot_name = name 
        for phenotype in phenotype_list:
            all_phenotype_list.append(phenotype)
            raw_interaction_list.append([prot_name, phenotype])
    # Delete phenotype if not present at least 2 times
    for interaction in raw_interaction_list:
        if all_phenotype_list.count(interaction[1]) >= 2:
            interaction_list.append(interaction)
    return interaction_list

def make_graph(interaction_list, filename):
    # Buidl graph
    G = nx.Graph()
    # Add all edge / nodes
    for interaction in interaction_list:
        G.add_edge(interaction[0], interaction[1])
    # Color
    color_map = []
    for node in G:
        if node in ["ADAMTS7", "ADAMTS12", "ADAMTS16", "ADAMTS18", "ADAMTS6", "ADAMTS10"]:
            color_map.append("blue")
        elif node in ["ADAMTS8", "ADAMTS1", "ADAMTS15", "ADAMTS5", "ADAMTS4", "ADAMTS9", "ADAMTS20"]:
            color_map.append("green")
        elif node in ["ADAMTS2", "ADAMTS3", "ADAMTS14"]:
            color_map.append("purple")
        elif node in ["ADAMTSL1", "ADAMTSL2", "ADAMTSL3", "ADAMTSL4", "ADAMTSL5", "ADAMTSL6/THDS4", "PAPLN", "ADAMTS13", "ADAMTS17", "ADAMTS19"]:
            color_map.append("red")
        else:
            color_map.append("grey")
    # Draw the graph
    nx.draw(G, with_labels=True, font_size=7, pos=nx.kamada_kawai_layout(G), node_color=color_map, edge_color="grey")
    plt.savefig(filename)
    plt.show()

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("leaf_file",
                        help = "csv file containing phenotype informations (presence of each phenotype for each leaf)",
                        type=str)
    parser.add_argument("equivalence_csv",
                        help = "csv file containing equivalence prot name and refseq",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    leaf_file = Path(args.leaf_file).resolve()
    equivalence_csv = args.equivalence_csv
    dict_name_refseq = load_equivalence(equivalence_csv)
    dict_prot_phenotypeList = parse_ppi_csv(leaf_file)
    interaction_list = make_interaction_list(dict_prot_phenotypeList, dict_name_refseq)
    make_graph(interaction_list, f"{leaf_file.parents[0]}/ppi_{leaf_file.stem}.pdf")


if __name__ == '__main__':
	main()

