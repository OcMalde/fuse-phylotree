#!/bin/python3

# Write proper csv leaf file from manual ppi file

import argparse
from pathlib import Path


def load_equivalence(equivalence_csv) -> dict:
    """
    Load equivalence csv file, 
    return a dict 
    { name : refseq }
    """
    dict_name_refseq = {}
    with open(equivalence_csv, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n","").split(",")
            name = splited_line[0]
            refseq = splited_line[1]
            dict_name_refseq[name] = refseq
    return dict_name_refseq

def load_ppi(ppi_csv) -> dict:
    """
    Load ppi csv file,
    csv format : 
    ppi_uniprot, ppi_name, my_prot_name
    return a dict
    { my_prot_name : [ppi list]}
    """
    dict_name_ppiList = {}
    with open(ppi_csv, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n","").split(",")
            name = splited_line[2]
            if len(splited_line[0]) > 1:
                ppi = splited_line[0]
            else:
                ppi = splited_line[1]
            if name in dict_name_ppiList:
                dict_name_ppiList[name].append(ppi)
            else:
                dict_name_ppiList[name] = [ppi]
    return dict_name_ppiList

def load_ppi_N2_format(ppi_csv) -> dict:
    dict_name_ppiList = {}
    with open(ppi_csv, "r") as csv_file:
        for line in csv_file:
            splited_line = line.replace("\n","").split(",")
            # Header line
            if "Symbol" in splited_line:
                header = splited_line
                for name in header[2:]:
                    dict_name_ppiList[name] = []
            # Data lines
            else:
                ppi_uniprot = splited_line[0]
                ppi_name = splited_line[1].replace(" ", "")
                for ppi_index, ppi_pres in enumerate(splited_line[2:], start=2):
                    if ppi_pres == "1":
                        prot_name = header[ppi_index]
                        dict_name_ppiList[prot_name].append(f"{ppi_uniprot}_{ppi_name}")
    return dict_name_ppiList
  
def clean_from_uniq_ppi(dict_name_ppiList) -> dict:
    """
    Keep only ppi if there is at least 2 proteins that shares this ppi
    """
    all_ppi_list = []
    for name, ppi_list in dict_name_ppiList.items():
        for ppi in ppi_list:
            all_ppi_list.append(ppi)
    uniq_func_list = list(set([f for f in all_ppi_list if all_ppi_list.count(f) > 1]))
    print("Phenotype(s) present at least 2 times :\t", " ".join(uniq_func_list), "\n")
    uniq_dict_name_ppiList = {name : [ppi for ppi in ppi_list if ppi in uniq_func_list] for name, ppi_list in dict_name_ppiList.items()}
    return uniq_dict_name_ppiList
              
def load_refseq_from_fasta(fasta_file) -> list:
    """
    Read fasta file and extract proper refseq list from it
    Return thema as a list
    """
    refseq_list = []
    with open(fasta_file, "r") as f_file:
        for line in f_file:
            if ">" in line:
                refseq = line.split("_")[0].replace(">","").replace("P","P_")
                refseq_list.append(refseq)
    return refseq_list

def write_leaf_csv(dict_name_ppiList, dict_name_refseq, refseq_interest, filename) -> None:
    """
    Use equivalence and ppi to write leaf csv file
    write only ppi of proteins that are in the fasta of interest
    csv format ; 
    refseq, ppi1|ppi2|ppi3
    """
    with open(filename, "w+") as csv_file:
        for name, ppiList in dict_name_ppiList.items():
            refseq = dict_name_refseq[name]
            if refseq in refseq_interest:
                print(name, refseq, "\t-> ", " ".join(ppiList))
                csv_file.write(f"{refseq},{'|'.join(ppiList)}\n")

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("ppi_csv",
                        help = "csv file containing ppi_uniprot, ppi_name, my_prot_name",
                        type=str)
    parser.add_argument("equivalence_csv",
                        help = "csv file containing my_prot_name, refseq",
                        type=str)
    parser.add_argument("corresponding_fasta",
                        help = "fasta file we are going to study",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    ppi_csv, equivalence_csv, fasta_file = args.ppi_csv, args.equivalence_csv, args.corresponding_fasta
    ppi_csv = Path(ppi_csv).resolve()
    equivalence_csv = Path(equivalence_csv).resolve()
    fasta_file = Path(fasta_file).resolve()
    dict_name_refseq = load_equivalence(equivalence_csv)
    #dict_name_ppiList = load_ppi(ppi_csv)
    dict_name_ppiList = load_ppi_N2_format(ppi_csv)
    # Keep only a ppi if at leat 2 times
    # dict_name_ppiList = clean_from_uniq_ppi(dict_name_ppiList)
    refseq_interest = load_refseq_from_fasta(fasta_file)
    filename = Path(f"{ppi_csv.parents[0]}/leaf_Manual_{fasta_file.stem}.csv").resolve()
    write_leaf_csv(dict_name_ppiList, dict_name_refseq, refseq_interest, filename)



if __name__ == '__main__':
	main()

