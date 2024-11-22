
#!/bin/python3

import argparse
import os
from pathlib import Path
from itertools import groupby
import urllib.request, urllib.parse, urllib.error

#==============================================================================
# Read files
#==============================================================================

def fasta_fast_parse(fasta_file) -> list:
    """
    Fast fasta parser
    https://www.biostars.org/p/710/ modified
    """
    fh = open(fasta_file)
    seq_list = []
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        seq_list.append((headerStr, seq))
    return seq_list

def get_refseq(dict_refseq_fasta, seq_list) -> tuple:
    """
    Extract only refseq from seq_list
    where a seq is a tuple (header, sequence)
    """
    refseq_list = []
    for seq in seq_list:
        refseq = seq[0].split(" ")[0].split(".")[0]
        dict_refseq_fasta[refseq] = seq
        refseq_list.append(refseq)
    return dict_refseq_fasta, refseq_list

def makeAssocDict(assocF) -> dict:
    """
    Open an association file and return his dictionnary
    association file : taxid,assoc
    dictionnary : {taxid : assoc}
    """
    dic_taxid_assoc = {}
    with open(assocF, "r") as a_file:
        for line in a_file:
            dic_taxid_assoc[line.split(",")[0]] = line.split(",")[1].replace("\n","").split(".")[0]
    return dic_taxid_assoc

#==============================================================================
# Search refseq in gff database
#==============================================================================

def search_locus_gff(refseq_list, gff_file) -> dict:
    """
    Search gff file for our refseq of interest
    """
    dict_refseq_locus = {refseq : [] for refseq in refseq_list} 
    with open(gff_file) as gff_f:
        for line in gff_f:
            if line.startswith("#!"):
                continue
            elif line.startswith("##"):
                continue
            else:
                for refseq in refseq_list:
                    if refseq in line:
                        splited_line = line.strip().split("\t")
                        assert len(splited_line) == 9
                        locus_info = {}
                        locus_info["seqid"] = urllib.parse.unquote(splited_line[0])
                        locus_info["strand"] = urllib.parse.unquote(splited_line[6])
                        locus_info["start"] = urllib.parse.unquote(splited_line[3])
                        locus_info["end"] = urllib.parse.unquote(splited_line[4])
                        dict_refseq_locus[refseq].append(locus_info)
    return dict_refseq_locus

#==============================================================================
# Assign proteins locus
#==============================================================================

def make_dic_locus_protList(dict_refseq_locus) -> dict:
    """
    Compare the gff feature of all refseq,
    And build a dictionnary
    {locus : [proteins list]}
    """
    dic_locus_protList = {}
    loc_id = 0
    for a_refseq, a_locus_list in dict_refseq_locus.items():
        assigned = False
        for a_locus_info in a_locus_list:
            for b_locus_id, b_prot_list in dic_locus_protList.items():
                for b_locus_list in [p[1] for p in b_prot_list]:
                    for b_locus_info in b_locus_list:
                        if same_locus(a_locus_info, b_locus_info):
                            if a_refseq not in [p[0] for p  in b_prot_list]:
                                b_prot_list.append((a_refseq, a_locus_list))
                            assigned = True
                            break
        if assigned == False:
            if len(a_locus_list) == 0:
                print(f"No annootation found for {a_refseq} ...")
            else:
                dic_locus_protList[f"{a_locus_list[0]['seqid']}_{loc_id}"] = [(a_refseq, a_locus_list)]
                loc_id += 1
    return dic_locus_protList

def same_locus(locus_1, locus_2) -> bool:
    """
    Compare if 2 locus info dict overlapp
    That means they are on the same sequence, same strand and overlapp on their positions
    """
    same = False
    if locus_1["seqid"] == locus_2["seqid"]:
        if locus_1["strand"] == locus_2["strand"]:
            if locus_1["end"] >= locus_2["start"] and locus_2["end"] >= locus_1["start"]:
                same = True
    return same

#==============================================================================
# Write output
#==============================================================================

def write_fasta_directory(dict_taxid_dictlocus, dict_refseq_fasta, directory) -> None:
    """
    Write directory containing fasta files
    """
    all_para_file = open(f"{directory}/longest_isoform.fasta", "w+")
    csv_str = ""
    for taxid, dict_locus in dict_taxid_dictlocus.items():
        for locus, prot_list in dict_locus.items():
            locus_file = open(f"{directory}/{taxid}_{locus}.fasta", "w+")
            longest_iso = ("None", "")
            csv_str += f"{taxid},{locus},{','.join(p[0] for p in prot_list)}\n"
            for prot in prot_list:
                header, sequence = dict_refseq_fasta[prot[0]]
                # Name for reconciliation refseq no _ taxid
                header = f"{header.split(' ')[0].replace('P_','P')}_{taxid}"
                locus_file.write(f">{header}\n{sequence}\n")
                if len(sequence) > len(longest_iso[1]):
                    longest_iso = (header, sequence)
            all_para_file.write(f">{longest_iso[0]}\n{longest_iso[1]}\n")
        locus_file.close()
    all_para_file.close()
    with open("locus_nbProt.csv", "w") as o_file:
        o_file.write(csv_str)

def getFastaNCBI(ncbi_id) -> str:
    """
    Search with ncbi api, the sequence on fasta format of the corresponding ncbi id
    """
    entrez = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&&id="
    rettype = "&rettype=fasta"
    fasta = requests.get(url = (entrez + ncbi_id + rettype)).text
    return fasta

def writeSpParaFile(dict_taxid_dictlocus, fn) -> None:
    """
    Write a csv file with the number of paralogs of each species
    """
    with open(fn, "w") as o_file:
        for sp, dict_locus in dict_taxid_dictlocus.items():
            o_file.write(f"{sp},{len(dict_locus)}\n")

def write_itol_nb(dict_taxid_dictlocus, dict_taxid_refseq, dn) -> None:
    """
    Write itol barchart file(s) 
    with number of sequences (og) per species, and number of uniq locus (ie paralogs)
    """
    # Init
    seq_itol_string = f"DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Nb_sequences\nCOLOR,#9370DB\nDATA\n"
    para_itol_string = f"DATASET_SIMPLEBAR\nSEPARATOR COMMA\nDATASET_LABEL,Nb_paralogs\nCOLOR,#00FF00\nDATA\n"
    # Data
    for taxid, dict_locus in dict_taxid_dictlocus.items():
        seq_itol_string += f"{taxid},{len(dict_taxid_refseq[taxid])}"
        para_itol_string += f"{taxid},{len(dict_locus)}\n"
    # Write file
    with open(f"{dn}/itol_nb_sequences.txt", "w+") as s_itol_file:
        s_itol_file.write(seq_itol_string)
    with open(f"{dn}/itol_nb_paralogs.txt", "w+") as p_itol_file:
        p_itol_file.write(para_itol_string)
                          
#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_directory",
                        help = "Directory containing proteins fasta selection files, name in the format selec_taxid.fasta",
                        type=str)
    parser.add_argument("--assoc_file",
                        help = ".csv association file (taxid, db_name)",
                        type=str)
    parser.add_argument("--gff_directory",
                        help = "Directory containing genome at the gff format (filename must contains db_name of the assoc_file)",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    
    # Get parsers input
    args = parser()
    fasta_directory = Path(args.fasta_directory).resolve()
    assoc_taxid_dbnt = Path(args.assoc_file).resolve()
    gff_directory = Path(args.gff_directory).resolve()
    
    dic_taxid_specie = makeAssocDict(assoc_taxid_dbnt)

    dict_refseq_fasta = {}
    dict_taxid_refseq = {}
    dict_taxid_dictlocus = {}
    # Iterate on fasta directory (1 fasta per specie, 1 gff per specie)
    for fasta_file in os.listdir(fasta_directory):
        
        # Read fasta file
        fn_fasta = Path(f"{fasta_directory}/{fasta_file}").resolve()
        if fn_fasta.suffix in (".fasta", ".fa", ".fas"):
            seq_list  = fasta_fast_parse(fn_fasta)
            dict_refseq_fasta, refseq_list = get_refseq(dict_refseq_fasta, seq_list)
        else:
            continue

        # Get the associated genome gff file and create / load database
        taxid = fn_fasta.stem.split("_")[1]
        print(f"Begin with {taxid} ...")
        genome_name = dic_taxid_specie[taxid]
        for gff_file in os.listdir(gff_directory):
            fn_gff = Path(f"{gff_directory}/{gff_file}").resolve()
            if genome_name in fn_gff.name:
                # Search locus of our refseq in the gff
                print(f"Mapping Refseq with gff features ...")
                # Search locus directly from gff
                dict_refseq_locus = search_locus_gff(refseq_list, fn_gff)

                # Compare features of genes and regroup per common locus
                print(f"Regroup Refseq per common locus ...")
                dict_locus_protList = make_dic_locus_protList(dict_refseq_locus)
                dict_taxid_dictlocus[taxid] = dict_locus_protList
                dict_taxid_refseq[taxid] = refseq_list

                break

    directory = Path(f"{fasta_directory}/isoforms_per_locus").resolve()
    if not os.path.exists(directory):
        os.makedirs(directory)
    writeSpParaFile(dict_taxid_dictlocus, f"{directory}/taxid_locusNB_{fasta_directory.name}.csv")
    write_itol_nb(dict_taxid_dictlocus, dict_taxid_refseq, directory)
    write_fasta_directory(dict_taxid_dictlocus, dict_refseq_fasta, directory)


if __name__ == "__main__":  
    main()

