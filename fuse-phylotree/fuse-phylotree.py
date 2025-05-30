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

    # Test for gene provide or not 
    if args.infer_gene_tree is None and args.gene_tree is None : 
        raise ValueError("You must provide a tree file if --infer_gene_tree is not used.")
    if args.infer_gene_tree is not None and args.gene_tree is not None : 
        raise ValueError("You use both the --infer_gene_tree option and supply a tree, choose one way only please")
    if args.leaf_functions_csv is not None and args.user_pastml_csv is not None : 
        print("You use both the --user_pastml_csv option and supply an annotation csv, the --user_pastml_csv will be used !")
    # Check if leaf_functions_csv looks like a tree file
    if args.leaf_functions_csv and args.leaf_functions_csv.endswith((".nwk", ".tree", ".newick")) and not args.gene_tree:
        # The user probably meant to skip leaf_functions_csv
        args.gene_tree = args.leaf_functions_csv
        args.leaf_functions_csv = None
    # Then validate:
    if not args.leaf_functions_csv and not args.user_pastml_csv:
        parser.error("You must provide either leaf_functions_csv or --user_pastml_csv.")

    # Prepare work directory
    if args.output_directory:
        work_dir = str(args.output_directory)
    else:
        work_dir = f"working_dir_{time.strftime('%Y%m%d-%H%M%S')}"
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
    if args.user_pastml_csv:
        user_pastml_csv = Path(args.user_pastml_csv).resolve()
        w_user_pastml_csv = Path(f"{work_dir}/{user_pastml_csv.name}").resolve()
        shutil.copy(user_pastml_csv, w_user_pastml_csv)
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

    # Gene tree
    if args.infer_gene_tree:
        # Launch the gene phylogeny
        print(f"[Gene tree inference - CAUTION ABOUT THE DEFAULT ROOT !!!] Begin gene phylogeny for {fasta_file.name} ...")
        gene_phylo_process, gene_tree_fn = gene_phylo.whole_phylo(w_fasta_file, species_tree)
        # Wait the gene phylo (need the gene tree to continue)
        gene_phylo_process.wait()
        print(f"Gene phylogeny for {fasta_file.name} is finished -> {gene_tree_fn.name}")
    else:
        # Use the provided gene tree
        gene_tree_fn = w_gene_tree
        print(f"Gene tree {gene_tree_fn.name} is provided ...")

    # Launch knwon domain / motifs recuperation
    print(f"Begin search of known domains / motifs for {fasta_file.name} ...")
    if args.reconc_domains:
        print(f"And domains fasta building")
        domains_segm_process_list, domains_fasta_tree_dn, domain_csv = get_domains.domains_fasta_aln_phylo(w_fasta_file)
        domain_process = subprocess.Popen("echo -n", shell=True, stdout=open(os.devnull, 'wb'))
    else:    
        domain_process, domain_csv = tools.known_domains(w_fasta_file)

    # Modules evolutions
    # Store module:gene mappings
    all_dict_gene_moduleList = []
    # Start module evolution iteration (we want to do it multiple times to have a more robust module:gene mapping)
    for i_module_evo in range(args.iter):

        # Launch the Modules segmentation
        extra_phyml_args = args.phyml_args.split() if args.phyml_args else None
        if args.plma_file:
            print(f"Paloma-D output {w_plma_output.name} is provided")
            modules_segm_process_list, modules_fasta_tree_dn = modules_segm.only_modules_phylo(w_fasta_file, w_plma_output, m_iter=i_module_evo, extra_args_phyml=extra_phyml_args)
        else:
            print(f"{i_module_evo+1}/{args.iter} Begin the modules segmentation for {fasta_file.name} ...")
            extra_paloma_args = args.paloma_args.split() if args.paloma_args else None
            modules_segm_process_list, modules_fasta_tree_dn = modules_segm.segmentation_and_modules_phylo(w_fasta_file, m_iter=i_module_evo, extra_args_paloma=extra_paloma_args, extra_args_phyml=extra_phyml_args)
        
        # Wait the modules segmentation (need the modules trees to continue)
        for ms_proc in modules_segm_process_list: 
            ms_proc.wait()
        if args.reconc_domains:
            [p.wait() for p in domains_segm_process_list]
            correct_domains_tree_process_list, correct_domains_tree_path_fn = get_domains.correct_domains_tree(domains_fasta_tree_dn, gene_tree_fn)
        print(f"{i_module_evo+1}/{args.iter} Modules segmentation for {fasta_file.name} is finished -> {modules_fasta_tree_dn.name}")
        # Write the main output with module descriptions
        modules_segm.write_2_module_descriptions(modules_fasta_tree_dn, f'{w_fasta_file.parents[1]}/2_module_descriptions.csv')
        print(f"{i_module_evo+1}/{args.iter} Begin modules tree correction using {gene_tree_fn.name} ...")
        extra_treefix_args = args.treefix_args.split() if args.treefix_args else None
        extra_raxml_args = args.raxml_args.split() if args.raxml_args else None
        correct_modules_tree_process_list, correct_modules_tree_path_fn = modules_segm.correct_modules_tree(modules_fasta_tree_dn, gene_tree_fn, extra_args_treefix=extra_treefix_args, extra_args_raxml=extra_raxml_args)
        # Wait for the end of these corrections (need it to continue)
        for tf_proc in correct_modules_tree_process_list: 
            tf_proc.wait()
        if args.reconc_domains:
            [p.wait() for p in correct_domains_tree_process_list]
            subprocess.run(f"cat {correct_domains_tree_path_fn} >> {correct_modules_tree_path_fn}", shell=True)
        
        print(f"{i_module_evo+1}/{args.iter} Modules trees correction finished -> {correct_modules_tree_path_fn.name}")
        
        # Make the DGS-reconciliation
        print(f"{i_module_evo+1}/{args.iter} Begin DGS reconciliation ...")
        extra_seadog_args = args.seadog_args.split() if args.seadog_args else None
        dgs_reconciliation_process, dgs_reconciliation_output_fn = tools.seadog_md(species_tree, gene_tree_fn, correct_modules_tree_path_fn, m_iter=i_module_evo, extra_args=extra_seadog_args)
        # Wait for dgs reconciliation to end (need the gene-specie reconciliation)
        dgs_reconciliation_process.wait()
        print(f"{i_module_evo+1}/{args.iter} DGS reconciliation finished")
        # Get the species-gene and the gene tree reconciliation events
        reconc_gene_tree, sp_gene_event_csv = integrate_3phylo.write_sp_gene_event(dgs_reconciliation_output_fn, gene_tree_fn)
        # TODO: is reconcilied gene tree always the same ? is sp/gene event always the same (fixed gene/species trees)

        # We only want to go until module:gene mappings (modules list associated to every genes)
        seadog_output = Path(dgs_reconciliation_output_fn).resolve()
        gene_tree_file = Path(gene_tree_fn).resolve()
        gene_tree, dict_module_mappingList, dict_module_modTree, gs_mapping_list = integrate_3phylo.read_seadogO(seadog_output, gene_tree_file)
        dict_gene_moduleList, dict_module_mappingList = integrate_3phylo.infers_modulesCompo(gene_tree, dict_module_mappingList, dict_module_modTree)    
        # Store module:gene mappings
        all_dict_gene_moduleList.append(dict_gene_moduleList)
    
    # Give the i produced mapping list directly to regroup them into one (filtering out every gene:modules not seen enough)
    module_list_csv = integrate_3phylo.rgrp_all_gene_module_lists(all_dict_gene_moduleList, freq_thres=args.mf_thres)

    # Launch the Ancestral scenario inference
    print(f"Begin ancestral scenario inference for {leaf_functions_csv.name} and {gene_tree_fn.name} ...")
    extra_pastml_args = args.pastml_args.split() if args.pastml_args else None
    w_user_pastml_csv = w_user_pastml_csv if args.user_pastml_csv else None
    ances_scenario_process, pastML_tab_fn = ances_scenario.acs_inference(reconc_gene_tree, w_leaf_functions_csv, sp_gene_event_csv, extra_args=extra_pastml_args, user_pastml_csv=w_user_pastml_csv)
    # Wait for ancestral scenario
    ances_scenario_process.wait()
    print(f"Finished the ancestral scenario inference -> {pastML_tab_fn.name}")
    # Wait for known domain search
    domain_process.wait()
    
    # Integrate all
    print(f"Begin DGS / gene and species phylo / function ancestral scenario integration ...")
    integrate_cmd = f"python3 {os.path.dirname(os.path.abspath(__file__))}/integrate_3phylo.py {dgs_reconciliation_output_fn} {gene_tree_fn} --pastml_tab {pastML_tab_fn} --domains_csv {domain_csv} --module_compositions {module_list_csv} --itol"
    integrate_process = subprocess.Popen(integrate_cmd, shell=True, stdout=open(os.devnull, 'wb'))
    integrate_process.wait()
    os.chdir(cwd)
    print(f"Finished integration")

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()

    # Fuse-phylotree arguments
    parser.add_argument("multi_fasta_file",
                        help = "Multi fasta file, with specific formated header >RefSeq_taxid (ex : >XP_012810820.2_8364)",
                        type=str)
    parser.add_argument("leaf_functions_csv",
                        help = "csv file containing for each of our sequence, the list of his functions (ex : XP_012810820.2, P59509 | P999999)",
                        type=str, nargs = "?")
    parser.add_argument("gene_tree",
                        help = "Gene tree to use as a support for the pastML and DGS reconciliation inference (WARNING, must correspond to the sequences in the multi fasta file !)",
                        type=str, nargs = "?")
    parser.add_argument("--output_directory",
                        help = "output directory name",
                        type=str)
    parser.add_argument("--iter",
                        help = "Number of times the whole module evolution inference will be performed, ie: module tree inference; their corrections; DGS reconciliation (default: 10)",
                        type=int,
                        default=10)
    parser.add_argument("--mf_thres",
                        help = "Module frequency thresold: frequency needed to consider a module is actually present at an ancestral gene over the module evolution iterations (default: 0.5)",
                        type=float,
                        default=0.5)
    parser.add_argument("--species_tree",
                        help = "Species tree to use as a support for the reconciliations (WARNING, must correspond to the taxid use in the other files !)",
                        type=str)
    parser.add_argument("--infer_gene_tree",
                        help = "Infer gene tree to use as a support for the pastML and DGS reconciliation inference (WARNING, user should check it and reroot it - we advise to only use it if you know what you are doing !)",
                        type=str)
    parser.add_argument("--plma_file",
                        help = "Paloma-2 output file (.agraph format, .dot, or .oplma format)",
                        type=str)
    parser.add_argument("--user_pastml_csv",
                        help = "PastML full input file, corresponding full custom states to use for the different sequence id (.csv format); eg, header: 'id,P59509,P999999', data: 'XP_012810820.2,1,0' or 'NP_001278744.1,0,,' ; unknown states (empty) will be inferred based on known states; sequence id will be converted to fit the reconcilied gene tree ids; ",
                        type=str)
    parser.add_argument("--reconc_domains",
                        help = "Do a DGS reconciliation with known modules (pfam / prosite) ; not tested",
                        action="store_true")
    
    # Custom third party software arguments (as raw strings)
    parser.add_argument("--paloma_args", 
                        type=str, 
                        help='Custom arguments to pass to paloma-D (e.g, "--thr 5 --min-size 5")')
    parser.add_argument("--phyml_args", 
                        type=str, 
                        help='Custom arguments to pass to PhyML for module trees inference (e.g, "--model JTT")')
    parser.add_argument("--treefix_args", 
                        type=str, 
                        help='Custom arguments to pass to TreeFix for gene-modules (e.g, "--niter 100 -D 1 -L 1" - corresponds to options inside -E from treefix or to --niter)')
    parser.add_argument("--raxml_args", 
                        type=str, 
                        help='Custom arguments to pass to RaxML (for TreeFix) for gene-modules (e.g, "-m PROTGAMMAJTT" - corresponds at options inside -e from treefix)')
    parser.add_argument("--seadog_args", 
                        type=str, 
                        help='Custom arguments to pass to SEADOG-MD (e.g, "--DD 5 --DL 1 --DTA 20 --GD 2 --GL 1")')
    parser.add_argument("--pastml_args", 
                        type=str, 
                        help='Custom arguments to pass to PastML (e.g, "--prediction_method ACCTRAN -m JTT")')
    
    return parser.parse_args()

#==============================================================================
# Main 
#==============================================================================

if __name__ == "__main__":
    # Get argparse arguments
    args = parser()
    # Call the phylo_char_mod function
    phylo_char_mod(args)



