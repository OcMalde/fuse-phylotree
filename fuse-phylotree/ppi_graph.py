#!/bin/python3

# Make a network x graph of protein-protein interactions from a psicquic

import argparse
import networkx as nx
import matplotlib.pyplot as plt


def parse_ppi_csv(ppi_file) -> list:
    """
    Read and parse a ppi csv file
    """
    ppi_list = []
    with open(ppi_file, "r") as csv_file:
        for line in csv_file:
            ppi_dict = {}
            splited_line = line.replace("\n","").split(",")
            if "refSeq_A" in splited_line: continue
            ppi_dict["refseq_A"] = splited_line[0]
            ppi_dict["uniprot_A"] = splited_line[1]
            ppi_dict["uniprot_B"] = splited_line[2]
            ppi_dict["species"] = splited_line[3]
            ppi_dict["database"] = splited_line[4]
            ppi_dict["evidence"] = splited_line[5]
            ppi_dict["pubmed"] = splited_line[6]
            ppi_list.append(ppi_dict)
    return ppi_list

def make_graph(ppi_list):
    G = nx.Graph()
    for ppi in ppi_list:
        G.add_edge(ppi["uniprot_A"], ppi["uniprot_B"])
    interactor_degree_2 = list(set([ppi["uniprot_B"] for ppi in ppi_list if G.degree[ppi["uniprot_B"]] > 1]))
    print(" ".join(interactor_degree_2))
    print(len(interactor_degree_2))
    
    nx.draw(G, with_labels=True, font_size=6, pos=nx.nx_agraph.graphviz_layout(G))
    plt.show()




#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("ppi_file",
                        help = "csv file containing protein-protein interaction informations, like the database they come from, or the publications",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    ppi_file = args.ppi_file
    ppi_list = parse_ppi_csv(ppi_file)
    make_graph(ppi_list)


if __name__ == '__main__':
	main()

