#!/bin/python3

# A tool to characterize proteics functionnal modules, using phylogeny

import argparse
import subprocess
import os
import time
import shutil
from pathlib import Path

import species_phylo
import gene_phylo
import modules_segm
import ances_scenario
import integrate_3phylo
import tools
import get_domains

#==============================================================================
# Phylo char mod function that run the whole programs
#==============================================================================

def phylo_char_mod(args) -> None:
    """
    Function that does the whole pipeline, 
    1) Module segmentation
    2) Gene phylogeny
    3) Ancestral scenario
    And integrate the whole to make modules <=> functions associations
    """
    print(f"Begin ...")
    fasta_file, leaf_functions_csv = args.multi_fasta_file, args.leaf_functions_csv
    cwd = os.getcwd()
    fasta_file = Path(fasta_file).resolve()
    leaf_functions_csv = Path(leaf_functions_csv).resolve()
    
    # Prepare work directory
    if args.output_directory:
        work_dir = str(args.output_directory)
    else:
        work_dir = f"multi_phylo_correlation_{fasta_file.stem}_{time.strftime('%Y%m%d-%H%M%S')}"
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    w_fasta_file = Path(f"{work_dir}/{fasta_file.name.replace('_','')}").resolve()
    w_leaf_functions_csv = Path(f"{work_dir}/{leaf_functions_csv.name}").resolve()
    shutil.copy(fasta_file, w_fasta_file)
    shutil.copy(leaf_functions_csv, w_leaf_functions_csv)
    if args.species_tree:
        species_tree = Path(args.species_tree).resolve()
        w_species_tree = Path(f"{work_dir}/{species_tree.name}").resolve()
        shutil.copy(species_tree, w_species_tree)
    if args.gene_tree:
        gene_tree = Path(args.gene_tree).resolve()
        w_gene_tree = Path(f"{work_dir}/{w_fasta_file.stem}.tree").resolve()
        shutil.copy(gene_tree, w_gene_tree)
    if args.plma_file:
        plma_output = Path(args.plma_file).resolve()
        w_plma_output = Path(f"{work_dir}/{w_fasta_file.stem}_{'_'.join(plma_output.name.split('_')[-2:])}").resolve()
        shutil.copy(plma_output, w_plma_output)
    os.chdir(work_dir)
    
    # Species phylogeny
    if args.species_tree:
        # Use the provided species tree
        species_tree = w_species_tree
        print(f"Species tree {species_tree.name} is provided ...")
    else:
        # First, get the species tree from NCBI taxonomy
        taxid_list = species_phylo.taxid_from_fasta(w_fasta_file)
        print(f"Recuperation of the species phylogeny for {' '.join(taxid_list)} ...")
        species_tree = species_phylo.getNCBItaxo(taxid_list)

    # Launch the Modules segmentation
    if args.plma_file:
        print(f"Paloma-2 output {w_plma_output.name} is provided")
        modules_segm_process_list, modules_fasta_tree_dn = modules_segm.only_modules_phylo(w_fasta_file, w_plma_output)
    else:
        print(f"Begin the modules segmentation for {fasta_file.name} ...")
        modules_segm_process_list, modules_fasta_tree_dn = modules_segm.segmentation_and_modules_phylo(w_fasta_file)
    
    # Launch knwon domain / motifs recuperation
    print(f"Begin search of known domains / motifs for {fasta_file.name} ...")
    if args.reconc_domains:
        print(f"And domains fasta building")
        domains_segm_process_list, domains_fasta_tree_dn, domain_csv = get_domains.domains_fasta_aln_phylo(w_fasta_file)
        domain_process = subprocess.Popen("echo -n", shell=True, stdout=open(os.devnull, 'wb'))
    else:    
        domain_process, domain_csv = tools.known_domains(w_fasta_file)
    
    # Gene tree
    if args.gene_tree:
        # Use the provided gene tree
        gene_tree_fn = w_gene_tree
        print(f"Gene tree {gene_tree_fn.name} is provided ...")
    else:
        # Launch the gene phylogeny
        print(f"Begin gene phylogeny for {fasta_file.name} ...")
        gene_phylo_process, gene_tree_fn = gene_phylo.whole_phylo(w_fasta_file, species_tree)
        # Wait the gene phylo (need the gene tree to continue)
        gene_phylo_process.wait()
        print(f"Gene phylogeny for {fasta_file.name} is finished -> {gene_tree_fn.name}")

    # Wait the modules segmentation (need the modules trees to continue)
    for ms_proc in modules_segm_process_list: 
        ms_proc.wait()
    if args.reconc_domains:
        [p.wait() for p in domains_segm_process_list]
        correct_domains_tree_process_list, correct_domains_tree_path_fn = get_domains.correct_domains_tree(domains_fasta_tree_dn, gene_tree_fn)
    print(f"Modules segmentation for {fasta_file.name} is finished -> {modules_fasta_tree_dn.name}")
    print(f"Begin modules tree correction using {gene_tree_fn.name} ...")
    correct_modules_tree_process_list, correct_modules_tree_path_fn = modules_segm.correct_modules_tree(modules_fasta_tree_dn, gene_tree_fn)
    # Wait for the end of these corrections (need it to continue)
    for tf_proc in correct_modules_tree_process_list: 
        tf_proc.wait()
    if args.reconc_domains:
        [p.wait() for p in correct_domains_tree_process_list]
        subprocess.run(f"cat {correct_domains_tree_path_fn} >> {correct_modules_tree_path_fn}", shell=True)
        
    print(f"Modules trees correction finished -> {correct_modules_tree_path_fn.name}")
    
    # Make the DGS-reconciliation
    print(f"Begin DGS reconciliation ...")
    dgs_reconciliation_process, dgs_reconciliation_output_fn = tools.seadog_md(species_tree, gene_tree_fn, correct_modules_tree_path_fn)
    # Wait for dgs reconciliation to end (need the gene-specie reconciliation)
    dgs_reconciliation_process.wait()
    print(f"DGS reconciliation finished")
    # Get the species-gene and the gene tree reconciliation events
    reconc_gene_tree, sp_gene_event_csv = integrate_3phylo.write_sp_gene_event(dgs_reconciliation_output_fn, gene_tree_fn)
    
    # Launch the Ancestral scenario inference
    print(f"Begin ancestral scenario inference for {leaf_functions_csv.name} and {gene_tree_fn.name} ...")
    ances_scenario_process, pastML_tab_fn = ances_scenario.acs_inference(reconc_gene_tree, w_leaf_functions_csv, sp_gene_event_csv)
    # Wait for ancestral scenario
    ances_scenario_process.wait()
    print(f"Finished the ancestral scenario inference -> {pastML_tab_fn.name}")
    # Wait for known domain search
    domain_process.wait()
    
    # Integrate all
    print(f"Begin DGS / gene and species phylo / function ancestral scenario integration ...")
    integrate_cmd = f"python3 {os.path.dirname(os.path.abspath(__file__))}/integrate_3phylo.py {dgs_reconciliation_output_fn} {gene_tree_fn} --pastml_tab {pastML_tab_fn} --domains_csv {domain_csv} --itol"
    integrate_process = subprocess.Popen(integrate_cmd, shell=True, stdout=open(os.devnull, 'wb'))
    integrate_process.wait()
    os.chdir(cwd)
    print(f"Finished integration")

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("multi_fasta_file",
                        help = "Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)",
                        type=str)
    parser.add_argument("leaf_functions_csv",
                        help = "csv file containing for each of our sequence, the list of his functions (ex : XP_012810820.2, P59509 | P999999)",
                        type=str)
    parser.add_argument("--output_directory",
                        help = "output directory name",
                        type=str)
    parser.add_argument("--species_tree",
                        help = "Species tree to use as a support for the reconciliations (WARNING, must correspond to the taxid use in the other files !)",
                        type=str)
    parser.add_argument("--gene_tree",
                        help = "Gene tree to use as a support for the pastML and DGS reconciliation inference (WARNING, must correspond to the sequences in the multi fasta file !)",
                        type=str)
    parser.add_argument("--plma_file",
                        help = "Paloma-2 output file (.agraph format, .dot, or .oplma format)",
                        type=str)
    parser.add_argument("--reconc_domains",
                        help = "Do a DGS reconciliation with known modules (pfam / prosite)",
                        action="store_true")
    return parser.parse_args()

#==============================================================================
# Main 
#==============================================================================

if __name__ == "__main__":
    # Get argparse arguments
    args = parser()
    # Call the phylo_char_mod function
    phylo_char_mod(args)



