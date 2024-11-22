# Build itol file for a expression file

import argparse
from pathlib import Path
from ete3 import Tree

#==============================================================================
# Read and load files
#==============================================================================

def binary_itol_dataset(csv_fn, tree_fn) -> dict:
    """
    Read the csv file containing the expression informations
    And the gene tree
    and write a binary itol file corresponding them
    """
    # Get leaves names
    tree = Tree(str(tree_fn), format=1)
    leaves_list = [node.name for node in tree.get_leaves()]
    # Get expression datas
    dict_leaf_data = {}
    with open(csv_fn, "r") as csv_file:
        # Iterate on lines
        for line in csv_file:
            splited_line = line.replace("\n","").split(",")
            # Header case
            if splited_line[0] == "protein":
                header = splited_line[2:]
            # Data(s) case(s)
            else:
                refseq = splited_line[1].replace("P_", "P")
                leaf = [l for l in leaves_list if refseq in l][0]
                dict_leaf_data[leaf] = splited_line[2:]
    # Build itol string
    itol_str = "DATASET_BINARY\nSEPARATOR COMMA\nDATASET_LABEL,proteins expressions\nCOLOR,#FF8C00\n"
    itol_str += f"FIELD_SHAPES,{','.join(['1']*len(header))}\nFIELD_LABELS,{','.join(header)}\n"
    itol_str += "DATA\n"
    # Data time
    for leaf, expression in dict_leaf_data.items():
        itol_str += f"{leaf},{','.join(expression)}\n"
    # Write itol binary file
    itol_fn = f"itol_{csv_fn.stem}.txt"
    with open(itol_fn, "w+") as itol_file:
        itol_file.write(itol_str)

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv_file",
                        help = "csv file containing the expression of our gene of interest (if possible with their refseq)",
                        type=str)
    parser.add_argument("--gene_tree",
                        help = "newick gene tree for the one we want to build an itol visualisation file",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    # Get arguments
    args = parser()
    csv_fn = Path(args.csv_file).resolve()
    tree_fn = Path(args.gene_tree).resolve()
    binary_itol_dataset(csv_fn, tree_fn)
    
    
if __name__ == "__main__":
    main()
