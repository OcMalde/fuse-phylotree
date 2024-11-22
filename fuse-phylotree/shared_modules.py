# Look for modules shared by proteins in a paloma decomposition

import os
import argparse
from pathlib import Path
from natsort import natsorted
from ete3 import Tree

#==============================================================================
# Read files
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

def get_fasta_modules(module_dir) -> dict:
    """
    Read all fasta of the module fastas directory
    return dict {module : {mod_header:seq}}
    mod_header contains module name, start ends and proteins name plus data name
    e.g. B692|1248|1254|_NP001361725.1_214
    """
    dict_module_protSeq = {}
    for file in Path(module_dir).iterdir():
        if file.suffix == ".fasta":
            dict_name_seq = get_fasta_from_file(file)
            dict_module_protSeq[file.stem] = dict_name_seq
    return dict_module_protSeq

#==============================================================================
# Search for modules shared by proteins of interest
#==============================================================================

def shared_mod(prot_list, dict_module_protSeq) -> tuple:
    """
    Return modules shared by all prot in prot list, 
    based on modules dict
    Also return list of proteins containing at least 1 of these modules
    """
    module_list = []
    one_mod_prot_list = []
    dict_prot_moduleList = {}
    for module, mod_dict in dict_module_protSeq.items():    
        seen = 0
        p_list = []
        for m_p_name, seq in mod_dict.items():
            p_name = m_p_name.split("_")[1]
            p_list.append(p_name)
            if p_name in prot_list:
                seen += 1
        if seen == len(prot_list):
            module_list.append(module)
            for p in p_list:
                if p not in one_mod_prot_list:
                    one_mod_prot_list.append(p)  
                    dict_prot_moduleList[p] = []
    for module, mod_dict in dict_module_protSeq.items(): 
        if module in module_list: 
            for m_p_name, seq in mod_dict.items():
                p_name = m_p_name.split("_")[1]
                m_name = module
                m_start = m_p_name.split("|")[1]
                m_end = m_p_name.split("|")[2]
                dict_prot_moduleList[p_name].append([m_name, m_start, m_end])
    return module_list, one_mod_prot_list, dict_prot_moduleList

#==============================================================================
# Write oplma signature
#==============================================================================

def write_oplma_modulesList(module_list, prot_list, dict_module_protSeq, dict_prot_seq, out_fn) -> None:
    """
    Write an oplma file for 1 module list (eg shared by proteins)
    Each module of this module list will be bloc in the oplma file
    """
    oplma_str = '<?xml version="1.0"?>\n<PLMA version="imitation v0 for module list">\n'
    oplma_str += '\t<resume>\n'
    oplma_str += f'\t\t<sequences size="{len(prot_list)}">\n'
    dict_desc_nb = {}
    nb = 1
    for prot in prot_list:
        oplma_str += f'\t\t\t<sequence tag="{prot}" data="{dict_prot_seq[prot]}" />\n'
        dict_desc_nb[prot] = nb
        nb += 1 
    oplma_str += '\t\t</sequences>\n'
    oplma_str += '\t</resume>\n'
    # Blocs/modules 
    oplma_str += '\t<partition>\n'
    for module in module_list:
        len_mod = len(list(dict_module_protSeq[module].values())[0])
        for res_index in range(len_mod):
            oplma_str += '\t\t<part>\n'
            for m_p_name in dict_module_protSeq[module]:
                p_name = m_p_name.split("_")[1]
                if p_name in prot_list:
                    m_start = int(m_p_name.split("|")[1])
                    pos = int(m_start) + int(res_index)
                    oplma_str += f'\t\t\t<aa seq="{dict_desc_nb[p_name]}" pos="{pos}" />\n' 
            oplma_str += '\t\t</part>\n'
    oplma_str += '\t</partition>\n'
    # End plma
    oplma_str += '</PLMA>\n'
    # Write file
    with open(out_fn, "w+") as oplma_file:
        oplma_file.write(oplma_str)

#==============================================================================
# Make protomata and logos visualisation pdf
#==============================================================================
    
def make_protomat_logo(oplma_fn) -> None:
    """
    Call protomata script to build protomat logo from  oplma file
    """
    os.system("rm -rf logo_dir")
    oplma_fn = Path(oplma_fn).resolve()
    os.system(f"~/Tools/protomata/bin/protobuild -i {oplma_fn} --proto-dot")
    proto_fn = f"{oplma_fn.stem}_q0.proto"
    os.system(f"~/Tools/protomata/scripts/proto2logo.py -i {proto_fn}")
    protodot_fn = f"{oplma_fn.stem}_q0_logos.dot"
    ps_fn = f"{oplma_fn.stem}_q0_logos.ps"
    os.system(f"dot -Tps2 logo_dir/{protodot_fn} -o {ps_fn}")
    os.system(f"ps2pdf {ps_fn}")
    #os.system(f"pstoedit -f plot-svg {ps_fn} {oplma_fn}_logos.svg")
    
def write_itol_shared(filename, prot_list, dict_prot_moduleList, dict_prot_seq, gene_tree, mod_list) -> None:
    with open(filename, "w+") as itol_file:
        itol_file.write("DATASET_DOMAINS\nSEPARATOR COMMA\nSHOW_INTERNAL,0\nSHOW_DOMAIN_LABELS,0\nBORDER_WIDTH,1\n")
        itol_file.write(f"DATASET_LABEL,{filename.split('.')[0]} shared\n")
        itol_file.write("WIDTH,1000\nHEIGHT_FACTOR,0.4\nBACKBONE_HEIGHT,0\n")
        itol_file.write("MARGIN,-1000\n")
        itol_file.write("COLOR,#ff0000\nDATA\n")   
        for gene in gene_tree.get_leaves():
            for prot, module_list in dict_prot_moduleList.items():
                max_pos = len(dict_prot_seq[prot])
                # If possess some modules (but not all)
                if prot in gene.name:
                    color = "#e5e500"
                    # if all modules cause one of the input
                    if prot in prot_list:
                        color = "#1E90FF"
                    # if all modules, but not an input
                    elif len(module_list) == len(mod_list):
                        color = "#FF7F50"
                    itol_file.write(f"{gene.name},{max_pos}")
                    for module in module_list:
                        m_name = module[0]
                        m_start = int(module[1]) - 1
                        m_end = module[2]
                        itol_file.write(f",{itol_shared(m_start, m_end, m_name, color)}")
                    itol_file.write("\n")

def itol_shared(start, end, name, color) -> str:
    return f"RE|{round(int(start),2)}|{round(int(end),2)}|{color}|{name}"

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("corresponding_fasta",
                        help = "fasta file with the full sequences of interest",
                        type=str)
    parser.add_argument("module_directory",
                        help = "directory containing fasta file for all our modules",
                        type=str)
    parser.add_argument("--protein_list_file",
                        help = "txt file containing a list of protein of interest (1 name per line), will search for the modules shared by all of them",
                        type=str)
    parser.add_argument("--gene_tree",
                        help = "gene tree file, permits to make itol file",
                        type=str)
    parser.add_argument("--at_least_1",
                        help = "show all proteins with at least 1 modules of the shared\nDefault being only proteins of interest",
                        action="store_true")
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    # Parsers
    args = parser()
    fasta_file = Path(args.corresponding_fasta).resolve()
    module_dir = Path(args.module_directory).resolve()
    if args.protein_list_file:
        txt_file = Path(args.protein_list_file)
        prot_list = []
        with open(txt_file, "r") as t_file:
            out_fn = f"shared_{txt_file.stem}.oplma"
            itol_fn = f"itolShared_{txt_file.stem}.txt"
            for line in t_file:
                prot_list.append(line.replace("\n",""))
    # Read and load all files
    dict_module_protSeq = get_fasta_modules(module_dir)
    dict_prot_seq_reconc = get_fasta_from_file(fasta_file)
    dict_prot_seq = {prot.split("_")[0] : seq for prot, seq in dict_prot_seq_reconc.items()}
    # Get shared modules
    print(f"Modules shared by {' '.join(prot_list)}:")
    module_list, one_mod_prot_list, dict_prot_moduleList = shared_mod(prot_list, dict_module_protSeq)
    module_list = natsorted(list(set(module_list)))
    print(f"{len(module_list)} modules\n{' '.join(module_list)}")
    # If we want all prot with at least 1 modules of the list (default being only modules on prot of interest)
    if args.at_least_1:
        prot_list = one_mod_prot_list
    # Write the oplma    
    write_oplma_modulesList(module_list, prot_list, dict_module_protSeq, dict_prot_seq, out_fn)
    make_protomat_logo(out_fn)
    if args.gene_tree:
        tree_file = Path(args.gene_tree).resolve()
        gene_tree = Tree(str(tree_file), format=1)
        write_itol_shared(itol_fn, prot_list, dict_prot_moduleList, dict_prot_seq, gene_tree, module_list)
 

if __name__ == '__main__':
	main()
    