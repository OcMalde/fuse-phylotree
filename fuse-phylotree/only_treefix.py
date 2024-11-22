#!/bin/python3


import os
import argparse
import subprocess
import shutil
import tools
from ete3 import Tree
from pathlib import Path

#==============================================================================
# Do the whole gene phylogeny
#==============================================================================

def treefix_tree(aln_file, species_tree, gene_tree, Verb=3, Niter=100) -> None:
    """
    Corrct my gene tree using species tree and his alignement file
    """
    # Treefix
    print(f"[Gene phylo] Begin tree correction ...")
    treefix_process, treefix_tree = tools.treefix(aln_file, gene_tree, species_tree, V=Verb, niter=Niter)
    # Treefix doesnt considere branch length, so we compute them on the fixed tree
    treefix_process.wait()
    print(f"[Gene phylo] Compute branch length of the fixed topology ...")
    compute_bl_process, treefix_bl_tree = tools.compute_branch_length(aln_file, treefix_tree)
    compute_bl_process.wait()
    # Reroot correctly the tree
    r_tree = ""
    with open(treefix_tree, "r") as t_file:
        for line in t_file:
            r_tree += line.replace("\n", "").replace(" ","")
    with open(treefix_tree, "w+") as t_file:
        t_file.write(r_tree)
    treefix_bl_tree = reroot_tree(treefix_tree, treefix_bl_tree) 

def reroot_tree(rooted_tree_fn, non_rooted_tree_fn) -> str:
    """
    Reroot a tree (with more info), based on a rooted tree with correct topology
    """
    rooted_tree = Tree(str(rooted_tree_fn), format=8)
    non_rooted_tree = Tree(str(non_rooted_tree_fn), format=2)
    non_rooted_tree.resolve_polytomy(recursive=True)
    new_tree = Path(f"{non_rooted_tree_fn.parents[0]}/final_{non_rooted_tree_fn.stem.split('.')[0].split('_')[-1]}.tree").resolve()
    non_rooted_tree.write(format=1, outfile=str(new_tree), format_root_node=True)
    return new_tree

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("aln_file",
                        help = "Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)",
                        type=str)
    parser.add_argument("gene_tree",
                        help = "Gene tree (newick format) to correct, corresponding the aln file",
                        type=str)
    parser.add_argument("species_tree",
                        help = "Species tree (newick format), of the taxid corresponding the sequences of the fasta file (ex : 8364)",
                        type=str)
    parser.add_argument("-V",
                        help = "Verbose level of treefix",
                        nargs='?',
                        const="3",
                        type=str)
    parser.add_argument("--niter",
                        help = "Number of iteration for treefix analyses",
                        nargs='?',
                        const="100",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    aln_file = Path(args.aln_file).resolve()
    species_tree = Path(args.species_tree).resolve()
    gene_tree = Path(args.gene_tree).resolve()

    treefix_tree(aln_file, species_tree, gene_tree, args.V, args.niter)


if __name__ == "__main__":
    main()

