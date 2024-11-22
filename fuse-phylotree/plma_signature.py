# Build a plma file for modules signatures (a "bloc" in the plma being a module in the signature)

import os
import argparse
from pathlib import Path
from ete3 import Tree

#==============================================================================
# Read files
#==============================================================================

def read_csv(csv_fn) -> dict:
    """ 
    Read csv file 
    format : gene,ppi_gain,ppi_lost,mod_gain,mod_lost
    return dict 
    {gene : ppi_gain, ppi_lost, mod_gain, mod_lost}
    """
    dict_gene_association = {}
    with open(csv_fn) as csv_file:
        for line in csv_file:
            if "gene" in line:
                continue
            else:
                splited_line = line.replace("\n","").split(",")
                gene = splited_line[0]
                ppi_gain_list = [i.split("_")[1] for i in splited_line[1].split("|") if len(i) > 1]
                ppi_lost_list = [i.split("_")[1] for i in splited_line[2].split("|") if len(i) > 1]
                mod_gain_list = [i for i in splited_line[3].split("|") if len(i) > 1]
                mod_lost_list = [i for i in splited_line[4].split("|") if len(i) > 1]
                dict_gene_association[gene] = [ppi_gain_list, ppi_lost_list, mod_gain_list, mod_lost_list]
    return dict_gene_association

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

def get_descendants(gene, gene_tree) -> list:
    """
    Return the descendants (only proteins name) of the gene, based on the gene tree topology
    """
    gene_node = gene_tree&gene
    desc_list = gene_node.get_leaves()
    desc_list = [desc.name.split("_")[0] for desc in desc_list]
    return desc_list

#==============================================================================
# Filters and operation
#==============================================================================

def co_appeared(dict_gene_association) -> dict:
    """
    Filter the gene association dict
    Keep only gene/association where at least 1 module and 1 ppi are gained
    Return subdict with same format
    """
    return {
        gene : association 
        for gene, association in dict_gene_association.items()
        if len(association[0]) > 0 and len(association[2]) > 0
        and gene.startswith("G")
        }

#==============================================================================
# Write oplma signature
#==============================================================================

def write_oplma_module(dict_module_protSeq, dict_prot_seq) -> None:
    """
    Write an oplma file for 1 module
    Each module will be the lonely bloc in the oplma file
    """
    out_dir = f"module_logo"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for module, prot_seq in dict_module_protSeq.items():
        oplma_str = '<?xml version="1.0"?>\n<PLMA version="imitation v0 for module">\n'
        oplma_str += '\t<resume>\n'
        len_mod = len(list(prot_seq.values())[0])
        oplma_str += f'\t\t<sequences size="{len_mod}">\n'
        dict_desc_nb = {}
        nb = 1
        for m_p_name, seq in prot_seq.items():
            p_name = m_p_name.split("_")[1]
            oplma_str += f'\t\t\t<sequence tag="{m_p_name}" data="{seq}" />\n'
            dict_desc_nb[m_p_name] = nb
            nb += 1
        oplma_str += '\t\t</sequences>\n'
        oplma_str += '\t</resume>\n'
        # Blocs/modules 
        oplma_str += '\t<partition>\n'
        for res_index in range(len_mod):
            oplma_str += '\t\t<part>\n'
            for m_p_name in prot_seq:
                pos = 1 + int(res_index)
                oplma_str += f'\t\t\t<aa seq="{dict_desc_nb[m_p_name]}" pos="{pos}" />\n' 
            oplma_str += '\t\t</part>\n'
        oplma_str += '\t</partition>\n'
        # End plma
        oplma_str += '</PLMA>\n'
        # Write file
        out_fn = Path(f"{out_dir}/{m_p_name.split('|')[0]}.oplma").resolve()
        with open(out_fn, "w+") as oplma_file:
            oplma_file.write(oplma_str)
        os.chdir(out_dir)
        make_protomat_logo(out_fn)
        os.chdir("..")
        
def write_oplma_SeqModules(seq_of_interest, gene_tree, dict_module_protSeq, dict_prot_seq, out_fn) -> None:
    """
    Write an oplma file for all modules of a selection of sequences
    Each module of this signature will be bloc in the oplma file
    """
    oplma_str = '<?xml version="1.0"?>\n<PLMA version="imitation v0 for module signatures">\n'
    oplma_str += '\t<resume>\n'
    # Get descendants and declare their sequences
    desc_list = [name for name, pos in seq_of_interest.items()]
    oplma_str += f'\t\t<sequences size="{len(desc_list)}">\n'
    dict_desc_nb = {}
    nb = 1
    for desc in desc_list:
        oplma_str += f'\t\t\t<sequence tag="{desc}" data="{dict_prot_seq[desc]}" />\n'
        dict_desc_nb[desc] = nb
        nb += 1
    oplma_str += '\t\t</sequences>\n'
    oplma_str += '\t</resume>\n'
    # Blocs/modules 
    oplma_str += '\t<partition>\n'
    module_list = []
    for module in dict_module_protSeq:
        len_mod = len(list(dict_module_protSeq[module].values())[0])
        for m_p_name in dict_module_protSeq[module]:
            p_name = m_p_name.split("_")[1]
            if p_name in desc_list:
                m_start = int(m_p_name.split("|")[1])
                if int(m_start) > int(seq_of_interest[p_name][0]) and int(m_start) < int(seq_of_interest[p_name][1]):
                    module_list.append(module)
    for module in module_list:
        len_mod = len(list(dict_module_protSeq[module].values())[0])
        #oplma_str += '\t\t<part>\n'
        for res_index in range(len_mod):
            oplma_str += '\t\t<part>\n'
            for m_p_name in dict_module_protSeq[module]:
                p_name = m_p_name.split("_")[1]
                if p_name in desc_list:
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
        
def write_dot_SeqModules(seq_of_interest, gene_tree, dict_module_protSeq, dict_prot_seq, out_fn) -> None:
    """
    Write a dot file for all modules of a selection of sequences
    Each module of this signature will be bloc in the dot file
    """
    # Declare graph and global parameters
    dot_str = 'DiGraph G\n{\n'
    dot_str += '\trankdir = LR;\n'
    dot_str += '\tconcentrate = true;\n'
    dot_str += '\tnodesep = 0.05;\n'
    dot_str += '\tranksep = 0.09;\n'
    dot_str += '\tnode [shape = record, fontname = monospace, height = 0.8, fontsize = 15];\n'
    dot_str += '\tedge [color = green, fontcolor = blue, weight = 1, headport = w, tailport = e, fontsize = 15];\n'
    subgraph_nb = 1
    # Get descendants and declare their sequences
    desc_list = [name for name, pos in seq_of_interest.items()]
    dot_str += f'\tsubgraph cluster_{subgraph_nb}\n'
    subgraph_nb += 1
    dot_str += "\t{\n"
    dot_str += "\t\tnode [shape = record, color = blue, fontcolor = black];\n"
    dot_str += "\t\tstyle = invis;\n"
    dict_desc_nb = {}
    nb = 1
    for desc in desc_list:
        dot_str += f'\t\t"({nb}, 1, 0)" [label = "{nb}: {desc}"];\n'
        dict_desc_nb[desc] = nb
        nb += 1
    dot_str += "\t}\n"
    # Declare end of the sequences
    dot_str += f'\tsubgraph cluster_{subgraph_nb}\n'
    subgraph_nb += 1
    dot_str += "\t{\n"
    dot_str += "\t\tnode [shape = record, color = green, fontcolor = green];\n"
    dot_str += "\t\tstyle = invis;\n"
    for desc in desc_list:
        nb = dict_desc_nb[desc]
        dot_str += f'\t\t"({nb}, {len(dict_prot_seq[desc])+1}, 0)" [label = "{nb}"];\n'
    dot_str += "\t}\n"   
    # Get Blocs/modules selection
    module_list = []
    dict_desc_moduleList = {desc : [] for desc in desc_list}
    for module in dict_module_protSeq:
        len_mod = len(list(dict_module_protSeq[module].values())[0])
        for m_p_name in dict_module_protSeq[module]:
            p_name = m_p_name.split("_")[1]
            if p_name in desc_list:
                m_start = int(m_p_name.split("|")[1])
                m_end = int(m_p_name.split("|")[2])
                if int(m_start) > int(seq_of_interest[p_name][0]) and int(m_start) < int(seq_of_interest[p_name][1]):
                    if module not in module_list:
                        module_list.append(module)
    # Declare blocs/modules subgraph
    for module in module_list:
        dot_str += f'\tsubgraph cluster_{subgraph_nb}\n'
        subgraph_nb += 1
        dot_str += "\t{\n"
        dot_str += "\t\tnode [shape = record, style = filled, color = grey, fontcolor = black];\n"
        dot_str += "\t\tcolor = chocolate1;\n"
        dot_str += f'\t\tlabel="{module}";\n'
        for m_p_name in dict_module_protSeq[module]:
            len_mod = len(list(dict_module_protSeq[module].values())[0])
            p_name = m_p_name.split("_")[1]
            if p_name in desc_list:
                nb = dict_desc_nb[p_name]
                m_seq = dict_module_protSeq[module][m_p_name]
                m_start = int(m_p_name.split("|")[1])
                m_end = int(m_p_name.split("|")[2])
                dot_str += f'\t\t"({nb}, {m_start}, {len(m_seq)})" [label = "{m_seq}"];\n'
                dict_desc_moduleList[p_name].append([m_start,m_end,len_mod])
        dot_str += "\t}\n"   
    # Declare all edges
    for desc in desc_list:
        m_list = dict_desc_moduleList[desc]
        m_list = sorted(m_list, key=lambda x: x[0])
        m_list.append([len(dict_prot_seq[desc])+1,len(dict_prot_seq[desc])+1,0])
        nb = dict_desc_nb[desc]
        prev_m = [1,1,0]
        for m in m_list:
            if m[0] == prev_m[1]+1:
                label = f'{nb}'
                color = "black"
                style = f'""'
            else:
                label = f'"{nb}: {prev_m[1]+1}..{m[0]-1}"'
                color = "black"
                style = f'""'
            dot_str += f'\t"({nb}, {prev_m[0]}, {prev_m[2]})" -> "({nb}, {m[0]}, {m[2]})" [label = {label}, color = {color}, style = {style}, fontcolor = blue, weight = 5];\n'	
            prev_m = m
    # End graph
    dot_str += "}\n"
    # Write file
    with open(out_fn, "w+") as dot_file:
        dot_file.write(dot_str)
        
def write_oplma_Allmodules(gene, gene_tree, dict_module_protSeq, dict_prot_seq, out_fn) -> None:
    """
    Write an oplma file for 1 module signature (ie module(s) gained at the given gene)
    Each module of this signature will be bloc in the oplma file
    """
    oplma_str = '<?xml version="1.0"?>\n<PLMA version="imitation v0 for module signatures">\n'
    oplma_str += '\t<resume>\n'
    # Get descendants and declare their sequences
    desc_list = get_descendants(gene, gene_tree)
    oplma_str += f'\t\t<sequences size="{len(desc_list)}">\n'
    dict_desc_nb = {}
    nb = 1
    for desc in desc_list:
        oplma_str += f'\t\t\t<sequence tag="{desc}" data="{dict_prot_seq[desc]}" />\n'
        dict_desc_nb[desc] = nb
        nb += 1
    oplma_str += '\t\t</sequences>\n'
    oplma_str += '\t</resume>\n'
    # Blocs/modules 
    oplma_str += '\t<partition>\n'
    for module in dict_module_protSeq:
        len_mod = len(list(dict_module_protSeq[module].values())[0])
        for res_index in range(len_mod):
            for m_p_name in dict_module_protSeq[module]:
                p_name = m_p_name.split("_")[1]
                if p_name in desc_list:
                    oplma_str += '\t\t<part>\n'
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

def write_oplma_signature(gene, module_list, gene_tree, dict_module_protSeq, dict_prot_seq, out_fn) -> None:
    """
    Write an oplma file for 1 module signature (ie module(s) gained at the given gene)
    Each module of this signature will be bloc in the oplma file
    """
    oplma_str = '<?xml version="1.0"?>\n<PLMA version="imitation v0 for module signatures">\n'
    oplma_str += '\t<resume>\n'
    # Get descendants and declare their sequences
    desc_list = get_descendants(gene, gene_tree)
    oplma_str += f'\t\t<sequences size="{len(desc_list)}">\n'
    dict_desc_nb = {}
    nb = 1
    for desc in desc_list:
        oplma_str += f'\t\t\t<sequence tag="{desc}" data="{dict_prot_seq[desc]}" />\n'
        dict_desc_nb[desc] = nb
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
                if p_name in desc_list:
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

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file",
                        help = "csv file containing the gene where modules and PPI co-appeared (format : gene,ppi_gain,ppi_lost,mod_gain,mod_lost, eg; G71_69_70,Q5XPI4_RNF123,,,B335|B1315)",
                        type=str)
    parser.add_argument("gene_tree",
                        help = "newick tree file containing the evolutionary relation shared by the gene in input csv_file",
                        type=str)
    parser.add_argument("corresponding_fasta",
                        help = "fasta file with the full sequences of interest",
                        type=str)
    parser.add_argument("module_directory",
                        help = "directory containing fasta file for all our modules",
                        type=str)
    parser.add_argument("--spec_refseq",
                        help = "csv file containing name of sequences of interest, with start and stop of region of interest, oplma of all their modules in this region will be produced",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    # Parsers
    args = parser()
    csv_file = Path(args.csv_file).resolve()
    tree_file = Path(args.gene_tree).resolve()
    fasta_file = Path(args.corresponding_fasta).resolve()
    module_dir = Path(args.module_directory).resolve()
    # Read and load all files
    dict_gene_association = read_csv(csv_file)
    gene_tree = Tree(str(tree_file), format=1)
    dict_module_protSeq = get_fasta_modules(module_dir)
    dict_prot_seq_reconc = get_fasta_from_file(fasta_file)
    dict_prot_seq = {prot.split("_")[0] : seq for prot, seq in dict_prot_seq_reconc.items()}
    # Keep only co-appeared at ancestral gene nodes
    dict_gene_association = co_appeared(dict_gene_association)
    # Write the oplma
    # For each module make a logo (long)
    #write_oplma_module(dict_module_protSeq, dict_prot_seq)
    
    if args.spec_refseq:
        refseq_file = Path(args.spec_refseq).resolve()
        seq_of_interest = {}
        with open(refseq_file) as t_file:
            for line in t_file:
                splited_line = line.replace("\n","").split(",")
                name = splited_line[0]
                start = splited_line[1]
                end = splited_line[2]
                seq_of_interest[name] = (start, end)
        write_dot_SeqModules(seq_of_interest, gene_tree, dict_module_protSeq, dict_prot_seq, f"specSeq_{refseq_file.stem}.dot")
        print(caca)
        write_oplma_SeqModules(seq_of_interest, gene_tree, dict_module_protSeq, dict_prot_seq, f"specSeq_{refseq_file.stem}.oplma")
        make_protomat_logo(f"specSeq_{refseq_file.stem}.oplma")
    
    for gene, association in dict_gene_association.items():
        module_list = association[2]
        write_oplma_signature(gene, module_list, gene_tree, dict_module_protSeq, dict_prot_seq, f"signature_{gene}.oplma")
        make_protomat_logo(f"signature_{gene}.oplma")
        #write_oplma_Allmodules(gene, gene_tree, dict_module_protSeq, dict_prot_seq, f"allModDesc_{gene}.oplma") 
        #make_protomat_logo(f"allModDesc_{gene}.oplma")
        # Human module signature (module only if present in at least 1 human sequence descendant)
        hs_module_list = []
        for module in module_list:
            for seq in dict_module_protSeq[module]:
                p_name = seq.split("_")[1]
                desc_list = get_descendants(gene, gene_tree)
                for prot, seq in dict_prot_seq_reconc.items():
                    if p_name in prot and p_name in desc_list:
                        if "9606" in prot and module not in hs_module_list:
                            hs_module_list.append(module)
                            continue
        print(gene, len(module_list), len(hs_module_list))
        write_oplma_signature(gene, hs_module_list, gene_tree, dict_module_protSeq, dict_prot_seq, f"hs_signature_{gene}.oplma")
        make_protomat_logo(f"hs_signature_{gene}.oplma")


if __name__ == '__main__':
	main()

