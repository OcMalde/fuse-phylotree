#!/bin/python3

# Infers an ancestral scenario for the known leafs functions

import os
import subprocess
import argparse
from ete3 import Tree
import tools
import shutil
from pathlib import Path

#==============================================================================
# Infers ancestral scenario from tree and function csv
#==============================================================================

def acs_inference(gene_tree_fn, gene_functions_csv, sp_gene_event_csv) -> tuple:
    """
    Ancestral scenario inference for a gene_tree and his leaf associated functions (from a csv)
    """
    # Prepare work directory
    acs_dir = f"{gene_tree_fn.parents[0]}/acs_dir_{gene_tree_fn.stem}"
    if not os.path.exists(acs_dir):
        os.makedirs(acs_dir)
    w_gene_tree_fn = Path(f"{acs_dir}/{gene_tree_fn.name}").resolve()
    w_gene_functions_csv = Path(f"{acs_dir}/{gene_functions_csv.name}").resolve()
    w_sp_gene_event_csv = Path(f"{acs_dir}/{sp_gene_event_csv.name}").resolve()
    shutil.copy(gene_tree_fn, w_gene_tree_fn)
    shutil.copy(gene_functions_csv, w_gene_functions_csv)
    shutil.copy(sp_gene_event_csv, w_sp_gene_event_csv)
    # Make the pastml input csv
    p_gene_tree_fn, pastml_csv = write_pastml_csv(w_gene_tree_fn, w_gene_functions_csv, w_sp_gene_event_csv)
    # Run pastml
    ances_scenario_process, pastML_tab_fn = tools.pastml(p_gene_tree_fn, pastml_csv)
    return ances_scenario_process, pastML_tab_fn

#==============================================================================
# Write the pastml input csv file
#==============================================================================

def write_pastml_csv(tree_file, leaf_functions_csv, sp_gene_event_csv) -> tuple:
    """
    Write a csv file that will be used as a pastml input
    https://github.com/evolbioinfo/pastml
    id, function1, function2
    leaf_node_name, (0|1||), (0|1||), ...
    """
    # Load the tree
    tree = Tree(str(tree_file), format=1)
    # Make a dict { leaf : [functions list]} and a uniq function list (init empty list for each gene)
    dict_leaf_funcList = {node.name : [] for node in tree}
    uniq_func_list = []
    with open(leaf_functions_csv, "r") as csv_file:
        for line in csv_file:
            splited_line= line.replace("\n","").replace(" ","").split(",")
            gene = splited_line[0]
            for node in tree.traverse("levelorder"):
                if gene.replace("_","") in node.name:
                    leaf = node.name 
            if 'leaf' in locals():
                functions_list = [f for f in splited_line[1].split('|') if len(f) > 1]
                dict_leaf_funcList[leaf] = functions_list
                for func in functions_list:
                    uniq_func_list.append(func.replace(' ',''))
                del leaf
    # Get the species-gene reconciliation event
    dict_gene_spGeneEvent = {}
    with open(sp_gene_event_csv, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n", "").split(",")
            gene = splited_line[0]
            for node in tree.traverse("levelorder"):
                if gene in node.name:
                    gene = node.name
            event = splited_line[1]
            dict_gene_spGeneEvent[gene] = event
        dict_gene_spGeneEvent[tree.get_tree_root().name] = ""
    # Search of in orthologs (orthologs after the last duplication of a gene)
    inOrtho_list = search_in_orthologs(dict_gene_spGeneEvent, tree)
    # Keep only phenotype if present at least 2 times (else no sens to infers acr)
    # uniq_func_list = list(set([f for f in uniq_func_list if uniq_func_list.count(f) > 1]))
    uniq_func_list = list(set([f for f in uniq_func_list if uniq_func_list.count(f) > 0]))
    # Make a presence / absence / dunno list for each leaf, stock all in dict with leaf keys 
    dict_leaf_pres = {}
    for leaf, func_list in dict_leaf_funcList.items():
        pres = []
        for func in uniq_func_list:
            if func in func_list:
                pres.append("1")
            else:
                # Check in orthologs, if they have the function, inferred it
                # ie : orthologs of a same gene duplication shared same function
                ortho_inferred = False
                # TEST
                """
                for inOrtho in inOrtho_list:
                    if leaf in inOrtho:
                        for ortho in inOrtho:
                            if func in dict_leaf_funcList[ortho] and ortho_inferred == False:
                                pres.append("1")
                                ortho_inferred = True
                """
                if ortho_inferred == False:
                    pres.append("0")
        # If the presence vector is empty, we consider we just dont have any informations about it
        if "1" not in pres:
            pres = [""] * len(pres) 
        assert len(pres) == len(uniq_func_list)
        # Add the binary profil of presence list (take in account the full 0000)
        dict_leaf_pres[leaf] = pres
        # Parse data ppi will could be by only adding gene with at least 1 1
        # So gene with any ppi will be inferred (false positiv risk)
    # Write patsml_csv
    pastml_csv = Path(f"{tree_file.parents[0]}/pastml_{tree_file.stem}_{leaf_functions_csv.stem}.csv").resolve()
    with open(pastml_csv, "w") as csv_file:
        csv_file.write(f"id,{','.join(uniq_func_list)}\n")
        for leaf, pres in dict_leaf_pres.items():
            csv_file.write(f"{leaf},{','.join(pres)}\n")
    return tree_file, pastml_csv

#==============================================================================
# Consider the species - genes reconciliations event
#==============================================================================

def search_in_orthologs(dict_gene_spGeneEvent, tree) -> list:
    """
    Search for in orthologs
    In orthologs are : gene common of a gene duplication events, only separate per spciation
    = all the genes after 'a last' gene duplication event
    = the same gene, only in different species = orthologs 1:1
    Return a list where each element is a list of these in orthologs genes
    (so each element represent a paralogs and his representants in diverses species)
    """
    # Init the list of the in ortho list
    inOrtho_list = []
    # Make list of genes that are gene duplications events
    geneDupli_list = [gene for gene, event in dict_gene_spGeneEvent.items() if event == "Gene duplication"]
    # Iterate on this genes (ancestral), to look at all their subtree
    for gene_duplication in geneDupli_list:
        # Get the subtree of this gene
        gene_duplication_node = tree&gene_duplication
        # Look at his parent, to see if its a speciation
        parent = gene_duplication_node.up
        if dict_gene_spGeneEvent[parent.name] == "Speciation":
            # Get the non gene_dupli child
            speciation_child = [child for child in parent.get_children() if dict_gene_spGeneEvent[child.name] in ("Speciation", "Leaf")]
            for sp_child in speciation_child:
                # Add his leaves in the dict (represent an paralogy 1:n)
                in_ortho = [leaf.name for leaf in sp_child.get_leaves()]
                if in_ortho not in inOrtho_list:
                    inOrtho_list.append(in_ortho)
        # Search his children (1 branch could have paralogy, the other not = n : 1)
        childs = gene_duplication_node.get_children()
        # Look at the 2 sides (2 childs, or more ...) of this gene duplication, 1 after another
        for child in childs:
            # Init the variables, a in para list, a out_para variable to known the situation
            out_para = False
            # Look at all their descendants
            for descendant in child.traverse("levelorder"):
                # Check his gene-species reconcilaition event
                # If its a gene duplication, were not in a in orthologs case
                if dict_gene_spGeneEvent[descendant.name] == "Gene duplication": 
                    out_para = True
            # If we"re not in an out_paralogy case
            if out_para == False:
                # We consider the leaf descendants of this gene duplication, as in orthologs
                inOrtho_list.append([leaf.name for leaf in child.get_leaves()])
    # Return the list of all the in orthologs list (1 per last gene duplication event)
    return inOrtho_list


