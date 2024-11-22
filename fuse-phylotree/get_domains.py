#!/bin/python3

import argparse
import os
import shutil
from Bio.ExPASy import ScanProsite
from prody import *
from pathlib import Path

import tools

confProDy(verbosity="critical")

#==============================================================================
# Domain Class
#==============================================================================
        
class protein:
    
    def __init__(self, header, sequence):
        self.header = header
        self.name = header.replace(">","")
        self.sequence = sequence
        self.domain_list = []
        self.addDomain("end", len(sequence), len(sequence))

    def addDomain(self, name, start, stop):
        domain = {}
        domain["name"] = name
        domain["start"] = start
        domain["stop"] = stop
        domain["sequence"] = self.sequence[start:stop+1] 
        self.domain_list.append(domain)

    def scanProsite(self):
        handle = ScanProsite.scan(seq=self.sequence, mirror='https://prosite.expasy.org')
        result = ScanProsite.read(handle)
        [self.addDomain(dom["signature_ac"], dom["start"], dom["stop"]) for dom in result]

    def scanDisorder(self):
        from iupred2a import iupred2a_lib
        result = iupred2a_lib.iupred(self.sequence.upper())
        [self.addDomain(f"Disorder_iupred2a_{result[res]}", res+1, res+1) for res in range(len(self.sequence.upper())) if result[res] > 0.5]

    def scanPfam(self):
        try :
            result = searchPfam(self.sequence)
            [self.addDomain(dom, int(result[dom]["locations"]["start"]), int(result[dom]["locations"]["end"])) 
                        for dom in result.keys()]
        except:
            print("Some fail on prody ...")

#==============================================================================
# Read multi fasta and load proteins
#==============================================================================

def get_fasta_from_file(multi_fasta) -> list:
    """
    Read the multi fasta file
    And parse it to return a protein list
    """
    protein_list = []
    sequence = ""
    with open(multi_fasta, "r") as fasta_file:
        for line in fasta_file:
            t_line = line.replace("\n","")
            if ">" in t_line:
                if sequence != "":
                    protein_list.append(protein(header, sequence))
                    sequence = ""
                header = t_line
            else: sequence += t_line
    protein_list.append(protein(header, sequence))
    return protein_list

def get_domains_fasta(protein_list, family_name) -> dict:
    """
    Regroup domains sequencces perr domains
    And  return them in a dict
    {domain : [sequence of this domainss list]}
    """
    dict_domain_domainSeqList = {}
    for protein in protein_list:
        protein_name = protein.name.split("_")[0]
        for domain in protein.domain_list:
            domain_name = domain["name"]
            start, end = domain["start"], domain["stop"]
            sequence = domain["sequence"]
            header = f"{domain_name}|{start}|{end}_{protein_name}_{family_name}"
            if len(sequence) > 1:
                if domain_name in dict_domain_domainSeqList:
                    dict_domain_domainSeqList[domain_name].append((header, sequence))
                else:
                    dict_domain_domainSeqList[domain_name] = [(header, sequence)]
    return dict_domain_domainSeqList

#==============================================================================
# Write csv output
#==============================================================================

def write_csv(protein_list, filename) -> None:
    """
    Write a csv output file
    protein, moduleName, start, stop
    """
    with open(filename, "w+") as csv_file:
        for prot in protein_list:
            for dom in prot.domain_list:
                csv_file.write(f"{prot.name},{dom['name']},{dom['start']},{dom['stop']}\n")
                
def write_domain_fastas(dict_domain_domainSeqList, family_name, directory) -> None:
    """
    Build a fasta file for each paloma bloc
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    for domain, seq_list in dict_domain_domainSeqList.items():
        if len(seq_list) == 1:
            continue
        with open(f"{directory}/{domain}.fasta","w") as f_file:
            for infos_seq in seq_list:
                header, seq = infos_seq[0], infos_seq[1]
                f_file.write(f">{header}\n{seq}\n")
        # If there is only 2 sequences in the blocs, make a trivial tree file
        if len(seq_list) == 2:
            names = []
            for infos_seq in seq_list:
                header = infos_seq[0]
                names.append(header)
            with open(f"{directory}/{domain}.tree","w") as t_file:
                t_file.write(f"({names[0]},{names[1]});")

#==============================================================================
# Read fasta, search for domains and then build domains fasta
#==============================================================================

def domains_fasta_aln_phylo(multi_fasta) -> tuple:
    """
    Read fasta, search for domains (pfam / proste)
    then build 1 fasta per domains (reconciliation format)
    """
    # Init files names
    multi_fasta = Path(multi_fasta).resolve()
    family_name = multi_fasta.stem.split("_")[0]
    out_fn = Path(f"{multi_fasta.parents[0]}/domains_{multi_fasta.stem}.csv").resolve()
    directory = Path(f"{multi_fasta.parents[0]}/domains_fasta").resolve()
    # Read fasta file
    protein_list = get_fasta_from_file(multi_fasta)
    # Prosite query
    [prot.scanProsite() for prot in protein_list]
    # Pfam query
    [prot.scanPfam() for prot in protein_list]
    # Write output csv file 
    write_csv(protein_list, out_fn)
    # Regroup per domains
    dict_domain_domainSeqList = get_domains_fasta(protein_list, family_name)
    # Write domains fasta files
    write_domain_fastas(dict_domain_domainSeqList, family_name, directory)
    # Go in the dir to launch alignement and phylo for each domains fasta file
    current = os.getcwd()
    os.chdir(directory)
    # Domains aligenements and phylogeny
    process_list = tools.all_msa_phylo(directory)
    os.chdir(current)
    # Return process, directory and csv filenames
    return process_list, directory, out_fn

#==============================================================================
# Correct domain tree, using gene tree
#==============================================================================

def correct_domains_tree(domains_fasta_tree_dn, gene_tree_fn) -> tuple:
    """
    Use treefix to correct domains trees of a directory, with treefix, using the gene tree
    """
    abs_path = ""
    process_list = []
    # Treefix for each domain
    for filename in Path(domains_fasta_tree_dn).iterdir():
        if filename.suffix == ".fasta" and filename.stem.startswith("muscle_"):
            fasta = Path(filename).resolve()
            tree = Path(f"{filename.parents[0]}/{filename.stem}.phylip_phyml_tree.txt").resolve()
            new_tree = Path(f"{filename.parents[0]}/{filename.stem}.tree").resolve()
            shutil.copy(tree, new_tree)
            treefix_process, treefix_tree = tools.treefix(fasta, new_tree, gene_tree_fn)
            process_list.append(treefix_process)
            abs_path += f"{os.path.abspath(treefix_tree)}\n"
    tree_path_fn = Path(f"domains_paths_{domains_fasta_tree_dn.stem}.txt").resolve()
    with open(tree_path_fn, "w") as text_file:
        text_file.write(abs_path)
    return process_list, tree_path_fn

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta",
                        help = "Multi fasta file, for the one we want to search de domains / motif presence",
                        type=str)
    parser.add_argument("-o",
                        help = "Output file name (.csv)",
                        type=str)
    parser.add_argument("--domain_fasta",
                        help = "Directory name where write domain fasta (1 fasta per domain)",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    # Get arguments
    args = parser()
    # Read fasta and load proteins
    multi_fasta = Path(args.fasta).resolve()
    protein_list = get_fasta_from_file(multi_fasta)
    # Prosite query
    [prot.scanProsite() for prot in protein_list]
    # Pfam query
    [prot.scanPfam() for prot in protein_list]
    # Disorder region prediction
    #[prot.scanDisorder() for prot in protein_list]
    # Write output csv file 
    write_csv(protein_list, args.o)
    
    if args.domain_fasta:
        family_name = multi_fasta.stem.split("_")[0]
        # Regroup per domains
        dict_domain_domainSeqList = get_domains_fasta(protein_list, family_name)
        # Write domains fasta files
        write_domain_fastas(dict_domain_domainSeqList, family_name, args.domain_fasta)


if __name__ == '__main__':
    main()
