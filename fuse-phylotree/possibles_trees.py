#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
from pathlib import Path
from ete3 import Tree
from itertools import cycle

#==============================================================================
# Load input
#==============================================================================

def read_trees(directory) -> list:
    """
    """
    tree_list = []
    for filename in directory.iterdir():
        fn = Path(filename).resolve()
        if fn.suffix in [".tree", ".nwk", ".nw", ".txt"]:
            tree = Tree(str(fn), format=1)
            tree_list.append(tree)
    return tree_list

def read_fasta(multi_fasta) -> dict:
    """
    """
    dict_leaf_fasta = {}
    sequence = ""
    with open(multi_fasta, "r") as fasta_file:
        for line in fasta_file:
            t_line = line.replace("\n","")
            if ">" in t_line:
                if sequence != "":
                    dict_leaf_fasta[leaf] = (header, sequence)
                    sequence = ""
                header = t_line
                leaf = t_line.replace(">","")
            else: sequence += t_line
    dict_leaf_fasta[leaf] = (header, sequence)
    return dict_leaf_fasta

#==============================================================================
# Combine trees
#==============================================================================

def combine_trees(tree_1, tree_2):
    """
    Combine 2 trees
        /--- Tree 1
    ---|
        \--- Tree 2
    """
    comb_tree = tree_1.copy("deepcopy")
    comb_tree.add_child(tree_2)
    comb_tree.set_outgroup(tree_2)
    return comb_tree

def build_all_possibles_trees(tree_list):
    """
    WARNING, ONLY FOR 3 TREES
    Concatenate all trees (all possibilities)
            /--- Tree 1
    ---|---|
       |    \--- Tree 2
       |
        \------- Tree 3
    And same with ((1,3),2); or ((2,3),1); 
    etc ...
    """
    tree_1, tree_2, tree_3 = tree_list
    # Case 1
    tree_12 = combine_trees(tree_1, tree_2)
    tree_123 = combine_trees(tree_12, tree_3)
    # Case 2
    tree_13 = combine_trees(tree_1, tree_3)
    tree_132 = combine_trees(tree_13, tree_2)
    # Case 3
    tree_23 = combine_trees(tree_2, tree_3)
    tree_231 = combine_trees(tree_23, tree_1)
    # Build list
    concat_trees = [tree_123, tree_132, tree_231, tree_12, tree_13, tree_23]
    # Return list with all combined trees 
    return concat_trees

#==============================================================================
# Write output
#==============================================================================

def build_trees_directory(tree_list, directory):
    """
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    i = 1
    for tree in tree_list:
        tree.write(format=1, outfile = f"{directory}/{i}_{len(tree.get_leaves())}_comb.tree", format_root_node=True)
        i += 1

def build_fasta(dict_leaf_fasta, tree_list, directory):
    """
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    i = 1
    for tree in tree_list:
        leaves = tree.get_leaves()
        with open(f"{directory}/{i}_{len(leaves)}_comb.fasta", "w") as fasta_file:
            for leaf in leaves:
                header, sequence = dict_leaf_fasta[leaf.name]
                fasta_file.write(f"{header}\n{sequence}\n")
        i += 1

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("newick_directory",
                        help = "Directory containing newick files of all the trees we want to combine",
                        type=str)
    parser.add_argument("--fasta_file",
                        help = "Multi fasta file, containing proteins sequences corresponding to our trees",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

if __name__ == "__main__":
    args = parser()
    directory = Path(args.newick_directory).resolve()
    tree_list = read_trees(directory)
    if len(tree_list) == 3:
        concat_trees = build_all_possibles_trees(tree_list) 
        out_directory = f"{directory.parents[0]}/combined_possibles_{directory.stem}"
        build_trees_directory(concat_trees, out_directory)
        if args.fasta_file:
            dict_leaf_fasta = read_fasta(Path(args.fasta_file).resolve())
            build_fasta(dict_leaf_fasta, concat_trees, out_directory)


