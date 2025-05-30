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

def segmentation_and_modules_phylo(fasta_file, m_iter, extra_args_paloma=None, extra_args_phyml=None) -> tuple:
    """
    Modules segmentation and the phylogeny for all modules
    """
    current = os.getcwd()
    # Prepare directory
    segm_dir = Path(f"{fasta_file.parents[0]}/modules_segm_dir_{fasta_file.stem}/").resolve()
    if not os.path.exists(segm_dir):
        os.makedirs(segm_dir)
    w_fasta_file = Path(f"{segm_dir}/{fasta_file.name}").resolve()
    shutil.copy(fasta_file, w_fasta_file)
    os.chdir(segm_dir)

    # Modules segmentation
    ms_process, ms_output = tools.segmentation(w_fasta_file, extra_args_paloma=extra_args_paloma, m_iter=m_iter)
    ms_process.wait()
    # Modules as fasta
    module_directory = tools.modules_fasta(ms_output)
    # Module evo iteration logic: 1 directory for the different trees
    m_iter_segm_dir = Path(f"{fasta_file.parents[0]}/modules_segm_dir_{fasta_file.stem}/{int(m_iter)+1}_mod_evo/").resolve()
    if not os.path.exists(m_iter_segm_dir):
        os.makedirs(m_iter_segm_dir)
    # Copy module_directory contents into m_iter_segm_dir
    for item in Path(module_directory).iterdir():
        dest = m_iter_segm_dir / item.name
        if item.is_dir():
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    # Change working directory to m_iter_segm_dir
    os.chdir(m_iter_segm_dir)
    # Modules phylogeny
    process_list = tools.all_phylo(m_iter_segm_dir, extra_args_phyml=extra_args_phyml)
    # Go back
    os.chdir(current)
    return process_list, m_iter_segm_dir

def only_modules_phylo(fasta_file, plma_output, m_iter, extra_args_phyml=None) -> tuple:
    """
    Modules phylogeny for already present module segmentation
    """
    current = os.getcwd()
    # Prepare directory
    segm_dir = Path(f"{fasta_file.parents[0]}/modules_segm_dir_{fasta_file.stem}/").resolve()
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
    # Module evo iteration logic: 1 directory for the different trees
    m_iter_segm_dir = Path(f"{fasta_file.parents[0]}/modules_segm_dir_{fasta_file.stem}/{int(m_iter)+1}_mod_evo/").resolve()
    if not os.path.exists(m_iter_segm_dir):
        os.makedirs(m_iter_segm_dir)
    # Copy module_directory contents into m_iter_segm_dir
    for item in Path(module_directory).iterdir():
        dest = m_iter_segm_dir / item.name
        if item.is_dir():
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    # Change working directory to m_iter_segm_dir
    os.chdir(m_iter_segm_dir)
    # Modules phylogeny
    process_list = tools.all_phylo(m_iter_segm_dir, extra_args_phyml=extra_args_phyml)
    os.chdir(current)
    return process_list, m_iter_segm_dir

#==============================================================================
# Correct modules tree, using gene tree
#==============================================================================

def correct_modules_tree(modules_fasta_tree_dn, gene_tree_fn, extra_args_treefix=None, extra_args_raxml=None) -> tuple:
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
            treefix_process, treefix_tree = tools.treefix(fasta, new_tree, gene_tree_fn, extra_args_treefix=extra_args_treefix, extra_args_raxml=extra_args_raxml)
            process_list.append(treefix_process)
            abs_path += f"{os.path.abspath(treefix_tree)}\n"
    tree_path_fn = Path(f"modules_paths_{modules_fasta_tree_dn.stem}.txt").resolve()
    with open(tree_path_fn, "w") as text_file:
        text_file.write(abs_path)
    return process_list, tree_path_fn

#==============================================================================
# Write main output with full module descriptions
#==============================================================================

def get_fasta_from_file(multi_fasta) -> dict:
    """
    Read the multi fasta file
    And parse it to return a protein dict
    """
    protein_dict = {}
    sequence = ""
    with open(multi_fasta, "r") as fasta_file:
        for line in fasta_file:
            t_line = line.replace("\n","")
            if ">" in t_line:
                if sequence != "":
                    protein_dict[header.replace(">","")] = sequence
                    sequence = ""
                header = t_line
            else: sequence += t_line.replace("X","A")
    protein_dict[header.replace(">","")] = sequence
    return protein_dict

def write_2_module_descriptions(modules_fasta_tree_dn, filename):
    """
    Write csv file with all module segments
    module_name, protein, start, end, segment
    """
    # Make the file 
    with open(filename, "w+") as csv_file:
        # Header
        csv_file.write("module,protein,start,end,segment\n") 
        # Look at all modules fasta (1 fasta = 1 module)
        for filename in Path(modules_fasta_tree_dn).iterdir():
            if filename.suffix == ".fasta":
                fasta = Path(filename).resolve()
                dict_name_seq = get_fasta_from_file(fasta)
                # For all segments in the module fasta file
                for name, segment in dict_name_seq.items():
                    module, start, end, node = name.split('|')
                    protein = node.split('_')[1]
                    # Write all line
                    csv_file.write(f"{module},{protein},{start},{end},{segment}\n")

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



