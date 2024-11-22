#!/bin/python3

import argparse
from ete3 import NCBITaxa, Tree, PhyloTree
from pathlib import Path

#==============================================================================
# Get NCBI taxonomy / build species tree
#==============================================================================

def getNCBItaxo(taxid_list) -> object:
    """
    Recuperation of the taxonomy topology of our species using the ncbi taxonomy database
    """
    ncbi = NCBITaxa()
    tree = ncbi.get_topology(taxid_list, intermediate_nodes=False)
    tree.resolve_polytomy(recursive=True)
    fn = Path(f"{len(taxid_list)}_species.tree").resolve()
    with open(fn, "w") as tree_file:
        tree_file.write(tree.write(format=8))
    return fn

def makeAssocDict(assocF) -> dict:
    """
    Open an association file and return his dictionnary
    association file : taxid,assoc
    dictionnary : {taxid : assoc}
    """
    dic_taxid_assoc = {}
    with open(assocF, "r") as a_file:
        for line in a_file:
            dic_taxid_assoc[line.split(",")[0]] = line.split(",")[1].replace("\n","")
    return dic_taxid_assoc

def taxid_from_fasta(fasta_file) -> list:
    """
    Get the taxid list from a fasta file
    where header are 
    >refseq_taxid
    """
    taxid_list = []
    with open(fasta_file, "r") as f_file:
        for line in f_file:
            if line.startswith(">"):
                taxid = line.replace("\n","").split("_")[1]
                if taxid not in taxid_list:
                    taxid_list.append(taxid)
    return taxid_list

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("taxid_file",
                        help = "Text file containing the list of the taxid (ncbi numbers) of the species where we want to search our orthologs",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    if args.taxid_file.endswith(".txt"):
        with open(args.taxid_file) as tax_file:
            taxid_list = [taxid for taxid in tax_file]
    elif args.taxid_file.endswith(".csv"):
        taxid_list = makeAssocDict(args.taxid_file).keys()
    getNCBItaxo(taxid_list)




if __name__ == '__main__':
	main()

