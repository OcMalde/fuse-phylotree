# Modules segmentation of a multi fasta file, and their corresponding phylogeny

import argparse
import os
import subprocess
import shutil
from pathlib import Path

import tools

#==============================================================================
# Segmentation and modules phylogeny
#==============================================================================

def segmentation_and_modules_phylo(fasta_file) -> tuple:
    """
    Modules segmentation and the phylogeny for all modules
    """
    current = os.getcwd()
    # Prepare directory
    segm_dir = Path(f"{fasta_file.parents[0]}/modules_segm_dir_{fasta_file.stem}").resolve()
    if not os.path.exists(segm_dir):
        os.makedirs(segm_dir)
    w_fasta_file = Path(f"{segm_dir}/{fasta_file.name}").resolve()
    shutil.copy(fasta_file, w_fasta_file)
    os.chdir(segm_dir)
    # Modules segmentation
    ms_process, ms_output = tools.segmentation(w_fasta_file)
    ms_process.wait()
    # Modules as fasta
    module_directory = tools.modules_fasta(ms_output)
    # Modules phylogeny
    process_list = tools.all_phylo(module_directory)
    os.chdir(current)
    return process_list, module_directory

def only_modules_phylo(fasta_file, plma_output) -> tuple:
    """
    Modules phylogeny for already present module segmentation
    """
    current = os.getcwd()
    # Prepare directory
    segm_dir = Path(f"{fasta_file.parents[0]}/modules_segm_dir_{fasta_file.stem}").resolve()
    if not os.path.exists(segm_dir):
        os.makedirs(segm_dir)   
    w_fasta_file = Path(f"{segm_dir}/{fasta_file.name}").resolve()
    shutil.copy(fasta_file, w_fasta_file)
    w_plma_output = Path(f"{segm_dir}/{plma_output.name}").resolve()
    shutil.copy(plma_output, w_plma_output)
    os.chdir(segm_dir)
    # Module segmentation is already done and provided
    ms_output = w_plma_output
    # Modules as fasta
    module_directory = tools.modules_fasta(ms_output)
    # Module phylogeny
    process_list = tools.all_phylo(module_directory)
    os.chdir(current)
    return process_list, module_directory

#==============================================================================
# Correct modules tree, using gene tree
#==============================================================================

def correct_modules_tree(modules_fasta_tree_dn, gene_tree_fn) -> tuple:
    """
    Use treefix to correct modules trees of a directory, with treefix, using the gene tree
    """
    abs_path = ""
    process_list = []
    # Treefix for each module
    for filename in Path(modules_fasta_tree_dn).iterdir():
        if filename.suffix == ".fasta":
            fasta = Path(filename).resolve()
            tree = Path(f"{filename.parents[0]}/{filename.stem}.phylip_phyml_tree.txt").resolve()
            new_tree = Path(f"{filename.parents[0]}/{filename.stem}.tree").resolve()
            shutil.copy(tree, new_tree)
            treefix_process, treefix_tree = tools.treefix(fasta, new_tree, gene_tree_fn)
            process_list.append(treefix_process)
            abs_path += f"{os.path.abspath(treefix_tree)}\n"
    tree_path_fn = Path(f"modules_paths_{modules_fasta_tree_dn.stem}.txt").resolve()
    with open(tree_path_fn, "w") as text_file:
        text_file.write(abs_path)
    return process_list, tree_path_fn

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("multi_fasta_file",
                        help = "Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()



