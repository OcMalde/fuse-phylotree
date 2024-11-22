#!/bin/python3

import argparse
from pathlib import Path


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
        }

def load_domains_from_csv(domains_file) -> dict:
    """
    Load the domains coordinates from a csv file
    csv file is : geneName (no node ID), domainsName, start, stop
    Return a dict
    { gene : [ domain_dict list ] }
    A domain of the domain_dict list is :
    { domain_name : [ start, end ]}
    """
    # Init the dict
    dict_gene_domainsList = {}
    # Open and read line per line the input csv_file
    with open(domains_file, "r") as csv_file:
        # Iterate on csv file
        for line in csv_file:
            # Split the line
            splited_line = line.replace("\n", "").split(",")
            g_name, d_name, d_start, d_end = splited_line[0], splited_line[1], splited_line[2], splited_line[3]
            # Get the gene name with node id, from the raw gene name of the csv file (nod node id)
            gene = g_name.split("_")[0]
            # Build the domains dict (1 csv line = 1 domain = 1 dict)
            dict_domain_posList = {f"{d_name}|{d_start}|{d_end}" : [d_start, d_end]}
            # Add it to the gene domains list dict
            if gene not in dict_gene_domainsList:
                dict_gene_domainsList[gene] = [dict_domain_posList]
            else:
                dict_gene_domainsList[gene].append(dict_domain_posList)
    # Return the dict
    return dict_gene_domainsList

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
    module header e.g. B692|1248|1254|_NP001361725.1_214
    return dict {module : {gene : [start, stop] } }
    """
    dict_module_genePos = {}
    for file in Path(module_dir).iterdir():
        if file.suffix == ".fasta":
            dict_name_seq = get_fasta_from_file(file)
            dict_gene_Mpos = {}
            for header in dict_name_seq:
                gene = header.split("_")[1]
                module = header.split("|")[0]
                start = header.split("|")[1]
                stop = header.split("|")[2]
                dict_gene_Mpos[gene] = [start, stop]
            dict_module_genePos[module] = dict_gene_Mpos
    return dict_module_genePos

def region_domains_TS(dict_gene_domainsList) -> dict:
    """
    Divide gene sequence by "regions", here based on domains
    {region : {gene : [start, stop] } }
    """
    dict_region_genePos = {}
    for gene, domains_list in dict_gene_domainsList.items():
        for domain_info in domains_list:
            dict_gene_Dpos = {}
            for d, info in domain_info.items():
                domain = d.split("|")[0]
                d_start = d.split("|")[1]
                d_stop = d.split("|")[2]
                dict_gene_Dpos[gene] = [d_start, d_stop]     
        dict_region_genePos[domain] = dict_gene_Dpos
    return dict_region_genePos

#==============================================================================
# Modules by region/domains
#==============================================================================

def modules_by_domains(module_list, dict_region_genePos, dict_module_genePos) -> dict:
    """
    Associate modules with region they are localised in 
    (If a modules, in 2 descendants is present in 2 different regions, we associate it with both)
    Return a dict 
    { region : [module list] }
    """
    dict_region_moduleList = {region : [] for region in dict_region_genePos}
    for region, gene_Rpos in dict_region_genePos.items():
        # By gene
        for gene, r_pos in gene_Rpos.items():
            r_start, r_stop = r_pos[0], r_pos[1]
            # Check for each module if its in this region
            for module in dict_module_genePos:
                m_info = dict_module_genePos[module]
                if gene in m_info:
                    m_pos = dict_module_genePos[module][gene]
                    m_start, m_stop = m_pos[0], m_pos[1]
                    # If in it, add it to region
                    if int(m_stop) >= int(r_start)-50 and int(r_stop)+50 >= int(m_start):
                        if module not in dict_region_moduleList[region]:
                            dict_region_moduleList[region].append(module)
    return dict_region_moduleList
                
#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file",
                        help = "csv file containing the gene where modules and PPI co-appeared (format : gene,ppi_gain,ppi_lost,mod_gain,mod_lost, eg; G71_69_70,Q5XPI4_RNF123,,,B335|B1315)",
                        type=str)
    parser.add_argument("--domains_csv",
                        help = "Add an domains csv file, at the format : gene_name, domain_name, start, end",
                        type=str)
    parser.add_argument("module_directory",
                        help = "directory containing fasta file for all our modules",
                        type=str)
    return parser.parse_args()



#==============================================================================
# Main
#==============================================================================

def main():
    # Get arguments
    args = parser()
    csv_fn = Path(args.csv_file).resolve()
    module_dir = Path(args.module_directory).resolve()
    dict_gene_association = read_csv(csv_fn)
    dict_gene_domainsList = load_domains_from_csv(args.domains_csv) 
    dict_module_genePos = get_fasta_modules(module_dir)
    dict_region_genePos = region_domains_TS(dict_gene_domainsList)
    
    
    print(f"{len(dict_gene_association)} module(s)-PPI(s) associations")
    dict_gene_association = co_appeared(dict_gene_association)
    print(f"{len(dict_gene_association)} co-apparitions associations")
    print(f"{len([gene for gene in dict_gene_association if gene.startswith('G')])} co-apparitions associations at ancestral gene node")
    


    # Count 
    ppi_list = []
    for gene, coapp in dict_gene_association.items():
        if gene.startswith("G"):
            for ppi in coapp[0]:
                if ppi not in ppi_list:
                    ppi_list.append(ppi)
    print(f"\t{len(ppi_list)} PPI associated with module signatures")
    mod_list = []
    for gene, coapp in dict_gene_association.items():
        if gene.startswith("G"):
            for mod in coapp[2]:
                if mod not in mod_list:
                    mod_list.append(mod)
    print(f"\t{len(mod_list)} modules present in these module signatures")
    
    print(mod_list)
    dict_region_moduleList = modules_by_domains(mod_list, dict_region_genePos, dict_module_genePos)
    for r, ml in dict_region_moduleList.items():
        print(r, len(ml), ml)
    
    # Select ancestral
    dict_gene_association = {gene : assoc for gene, assoc in dict_gene_association.items() if gene.startswith("G")}
    
    
    
    
if __name__ == "__main__":
    main()