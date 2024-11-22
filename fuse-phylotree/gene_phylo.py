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

def whole_phylo(fasta_file, species_tree) -> tuple:
    """
    Make the whole phylogeny for the given fasta file
    Using the corresponding species tree to correct it
    Return a tuple (process, gene_tree_fn)
    """
    # Prepare directory
    gene_dir = f"{fasta_file.parents[0]}/gene_phylo_dir_{fasta_file.stem}"
    if not os.path.exists(gene_dir):
        os.makedirs(gene_dir)
    w_fasta_file = Path(f"{gene_dir}/{fasta_file.name}").resolve()
    w_species_tree = Path(f"{gene_dir}/{species_tree.name}").resolve()
    shutil.copy(fasta_file, w_fasta_file)
    shutil.copy(species_tree, w_species_tree)
    # Alignement
    print(f"[Gene phylo] Begin MSA ...")
    msa_process, msa_fasta = tools.msa(w_fasta_file)
    msa_process.wait()
    print(f"[Gene phylo] MSA finished -> {msa_fasta.name}")
    # Trimal
    print(f"[Gene phylo] Begin MSA trimming ...")
    trimal_process, trimal_fasta = tools.trimal(msa_fasta)
    trimal_process.wait()
    print(f"[Gene phylo] MSA trimming finished -> {trimal_fasta.name}")
    # Phylogeny inference
    print(f"[Gene phylo] Begin phylogenetic inference ...")
    phylo_process, phylo_tree = tools.phylo(trimal_fasta)
    phylo_process.wait()
    clean_phylo_tree = Path(f"{trimal_fasta.parents[0]}/phyml_{trimal_fasta.stem}.tree").resolve()
    shutil.copy(phylo_tree, clean_phylo_tree)
    print(f"[Gene phylo] Phylogenetic inference finished -> {clean_phylo_tree.name}")
    # Treefix
    print(f"[Gene phylo] Begin tree correction ...")
    treefix_process, treefix_tree = tools.treefix(trimal_fasta, clean_phylo_tree, w_species_tree)
    # Treefix doesnt considere branch length, so we compute them on the fixed tree
    treefix_process.wait()
    print(f"[Gene phylo] Compute branch length of the fixed topology ...")
    compute_bl_process, treefix_bl_tree = tools.compute_branch_length(trimal_fasta, treefix_tree)
    compute_bl_process.wait()
    # Reroot correctly the tree
    r_tree = ""
    with open(treefix_tree, "r") as t_file:
        for line in t_file:
            r_tree += line.replace("\n", "").replace(" ","")
    with open(treefix_tree, "w+") as t_file:
        t_file.write(r_tree)
    treefix_bl_tree = reroot_tree(treefix_tree, treefix_bl_tree) 
    
    #treefix_bl_tree = treefix_tree

    # Return process
    return compute_bl_process, treefix_bl_tree

def reroot_tree(rooted_tree_fn, non_rooted_tree_fn) -> str:
    """
    Reroot a tree (with more info), based on a rooted tree with correct topology
    """
    rooted_tree = Tree(str(rooted_tree_fn), format=8)
    non_rooted_tree = Tree(str(non_rooted_tree_fn), format=2)
    non_rooted_tree.resolve_polytomy(recursive=True)
    new_tree = Path(f"{non_rooted_tree_fn.parents[0]}/{non_rooted_tree_fn.stem.split('.')[0].split('_')[-1]}.tree").resolve()
    non_rooted_tree.write(format=1, outfile=str(new_tree), format_root_node=True)
    return new_tree

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_file",
                        help = "Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)",
                        type=str)
    parser.add_argument("species_tree",
                        help = "Species tree (newick format), of the taxid corresponding the sequences of the fasta file (ex : 8364)",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    file_fasta = Path(args.fasta_file)
    tree_species = Path(args.species_tree)
    whole_phylo(file_fasta, tree_species)


if __name__ == "__main__":
    main()

