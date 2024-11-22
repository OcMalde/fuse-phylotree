#!/bin/python3

"""
Select homolog proteins from orthogroups containing a set of proteins of interest
"""

import os
import argparse
import re
import requests
from pathlib import Path
from tqdm import tqdm

#==============================================================================
# Read files
#==============================================================================

def read_orthogroups(fn) -> dict:
    """
    Read the Orthogroups csv file of orthofinder
    and return a dict
    { orthogroup_name : [ protein refseq list ]}
    
    Parameters
    ----------
    fn : str
        Name of a orthogroup file in csv format
        
    Returns
    -------
    dict_og_protList : dict
        Dictionary of orthogroups, orthgroup name (str) as key, list of refseq (str) clustered in this orthgroup as value
    """
    dict_og_protList = {}
    with open(fn, "r") as ortho_file:
        orthoGr_list = [ortho_group for ortho_group in ortho_file]
        for og in orthoGr_list:
            prot_in_group = re.split("\t|, ", og)
            if "Orthogroup" in prot_in_group:
                og_N = False
                continue
            elif "HOG" in prot_in_group: 
                og_N = True
                continue
            if og_N == False:
                og_name = prot_in_group[0]
                prot_list = [prot.replace("\n","") for prot in prot_in_group[1:] if prot not in ("", "\n")]
                dict_og_protList[og_name] = prot_list
            elif og_N == True:
                og_name = prot_in_group[1]
                prot_list = [prot.replace("\n","") for prot in prot_in_group[3:] if prot not in ("", "\n")]
                dict_og_protList[og_name] = prot_list
    return dict_og_protList

def read_prot_list(fn) -> list:
    """
    read a file with proteins of interest
    1 protein refseq /id per line
    
    Parameters
    ----------
    fn : str
        Name of a file in txt format (1 refseq id per line)
        
    Returns
    -------
    prot_list : list
        List of proteins refseq (str)
    """
    prot_list = []
    with open(fn, "r") as prot_file:
        for prot in prot_file : 
            if "." in prot:
                prot = prot.replace("\n", "")[:-2]
            else:
                prot = prot.replace("\n", "")
            prot_list.append(prot)   
    return prot_list

#==============================================================================
# Select my ortho / paralogs group
#==============================================================================

def select_my_family(dict_og_protList, prot_list) -> dict:
    """
    Select the orthogroups containing at least 1 of the protein of the prot list
    return a dict
    { orthogroup_name : [ protein refseq list ]}
    
    Parameters
    ----------
    dict_og_protList : dict
        Dictionry of orthogroups, orthgroup name (str) as key, list of refseq (str) clustered in this orthgroup as value
    prot_list : list
        List of proteins refseq (str)
        
    Returns
    -------
    my_og : dict
        Dictionary of selected orthogroups, orthgroup name (str) as key, list of refseq (str) clustered in this orthgroup as value
    """
    # Proteins in orthogroups
    # my_og = {
    #     og_name : og_prot_list 
    #     for og_name, og_prot_list in dict_og_protList.items() 
    #     if len(set(og_prot_list).intersection(set(prot_list))) > 0
    #     }
    my_og = {
        og_name: [prot[:-2] for prot in og_prot_list]
        for og_name, og_prot_list in dict_og_protList.items() 
        if any(prot[:-2] in prot_list for prot in og_prot_list)
        }
    # Proteins not in orthogroups
    seen_prot = []
    for og, p_list in my_og.items():
        for prot in p_list:
            if prot not in seen_prot: 
                seen_prot.append(prot)
    no_og = []
    for prot in prot_list:
        if prot not in seen_prot:
            no_og.append(prot)
    my_og["no_og"] = no_og
    return my_og

def makeAssocDict(assocF) -> dict:
    """
    Open an association file and return his dictionnary
    association file : taxid,assoc
    dictionnary : {taxid : assoc}
    
    Parameters
    ----------
    assocF : str
        Name of a association file in csv format (taxid,assoc)
        
    Returns
    -------
    dic_taxid_assoc : dict
        Dictionary of associations, taxid (str) as key, associated information (str) as value
    """
    dic_taxid_assoc = {}
    with open(assocF, "r") as a_file:
        for line in a_file:
            dic_taxid_assoc[line.split(",")[0]] = line.split(",")[1].replace("\n","")
    return dic_taxid_assoc

#==============================================================================
# Get fasta
#==============================================================================

def getFastaNCBI(ncbi_id) -> str:
    """
    Search with ncbi api, the sequence on fasta format of the corresponding ncbi id
    
    Parameters
    ----------
    ncbi_id : str
        Identifiant compatible with NCBI (e.g., refseq)
        
    Returns
    -------
    fasta : str
        Sequence in fasta format
    """
    entrez= "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&&id="
    rettype = "&rettype=fasta"
    fasta = requests.get(url = (entrez + ncbi_id + rettype)).text
    return fasta

#==============================================================================
# Write output
#==============================================================================

def write_orthogroups_csv(dict_og_protList, filename) -> None:
    """
    Write csv file, summering a orthogroups selection
    
    Parameters
    ----------
    dict_og_protList : dict
        Dictionary of orthogroups, orthgroup name (str) as key, list of refseq (str) clustered in this orthgroup as value
    filename : str
        Name to give to the output file, csv format
    """
    with open(filename, "w+") as csv_file:
        csv_file.write("og_name,refseq\n")
        for og, prot_list in dict_og_protList.items():
            for prot in prot_list:
                csv_file.write(f"{og},{prot}\n")
            

def writeFamilyFasta(my_family, dic_taxid_assoc, directory_name) -> None:
    """
    Write the multi fasta file with all the sequence of the proteins of the family
    
    Parameters
    ----------
    my_family : list
        List of selected proteins
    dic_taxid_assoc : dict
        Dictionary of associations, taxid (str) as key, associated information (str) as value
    directory_name : str
        Name of the output directory
    """
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    directory = Path(directory_name).resolve()
    print(f"Download the {len(my_family)} proteins of interest ...")
    for prot_id in tqdm(my_family):
        fasta_str = getFastaNCBI(prot_id)
        for taxid, specie in dic_taxid_assoc.items():
            filename = f"{directory}/selec_{taxid}.fasta"
            with open(filename, "a") as sp_fasta:
                if specie in fasta_str:
                    sp_fasta.write(fasta_str)

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("orthogroups_file",
                        help = "csv file containing the orthogroups from an orthofinder analysis",
                        type=str)
    parser.add_argument("myProtein_file",
                        help = "text file containing the ncbi id of the our protein of interest",
                        type=str)
    parser.add_argument("assocF_taxid_sp",
                        help = "association file : taxid, specieName",
                        type=str)
    parser.add_argument("--download",
                        help = "download proteins sequences of all the proteins in our orthogroups of interest",
                        action="store_true")
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    """
    Main function
    """
    # Get user arguments
    args = parser()
    og_file = Path(args.orthogroups_file).resolve()
    prot_file = Path(args.myProtein_file).resolve()
    assoc_file = Path(args.assocF_taxid_sp).resolve()
    # Read files 
    dict_og_protList = read_orthogroups(og_file)
    prot_list = read_prot_list(prot_file)
    dic_taxid_sp = makeAssocDict(assoc_file)
    # Select orthogroups with my protein of interest
    my_og = select_my_family(dict_og_protList, prot_list)
    # Some nice numbers
    print(f"*** My {len(prot_list)} proteins of interests are presents in {len(my_og)-1} orthogroups ***\n")
    for k, v in my_og.items():
        intersect = set(v).intersection(set(prot_list))
        print(f"=> {k} containing {len(v)} with {len(intersect)} of the proteins of interests")
        print(f"{' '.join(intersect)}\n")
    write_orthogroups_csv(my_og, f"selected_og_{prot_file.stem}.csv")
    if args.download:
        # Write fasta file with sequences of the selected orthogroups
        my_family_prot = []
        for k,v in my_og.items():
            for p in v: my_family_prot.append(p)
        out_directory = f"orthogroups_{prot_file.stem}"
        writeFamilyFasta(my_family_prot, dic_taxid_sp, out_directory)




if __name__ == '__main__':
	main()


