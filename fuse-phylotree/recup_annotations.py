#!/bin/python3

# Script for recuperation of protein-protein interaction in psicquic, for multiples proteins for 1 gene

import os
import pathlib
import time
import argparse
import requests
from urllib.parse import quote
from urllib.request import urlopen
import xml.etree.ElementTree as ET
from bioservices.uniprot import UniProt

#==============================================================================
# Read csv file
#==============================================================================

def read_locus_prots_csv(locus_protRefSeq_file) -> dict:
    """
    Read and load the information of the csv
    the csv must have the following format :
    taxid, genomic locus, refseq of proteins (separate by ,)
    ex : locus_nbProt.csv, build using spaln_prot_byLocus.py script
    return a dict at the format :
    {species : [list of gene]}
    where a gene is a dict :
    {gene : [list of refSeq]}
    """
    # Init the output dict
    dict_specie_geneList = {}
    # Read the csv file
    with open(locus_protRefSeq_file, "r") as csv_file:
        # Iterate on lines
        for line in csv_file:
            # Split it using the seperator, here ','
            splited_line = line.replace("\n","").split(",")
            specie = splited_line[0]
            locus = splited_line[1]
            refSeq_list = splited_line[2:]
            # Build the actual gene dict
            gene = {locus : refSeq_list}
            # Add it to the dict
            if specie not in dict_specie_geneList.keys():
                dict_specie_geneList[specie] = [gene]
            else:
                dict_specie_geneList[specie].append(gene)
    return dict_specie_geneList

#==============================================================================
# Prepare list to query
#==============================================================================

def make_query_list(dict_specie_geneList) -> list:
    """
    Build a list with tuple for query
    (refseq id, specie taxid)
    """
    refSeq_sp_list = []
    for sp, geneList in dict_specie_geneList.items():
        for gene in geneList:
            for locus, refSeq_list in gene.items():
                for refSeq in refSeq_list:
                    refSeq_sp_list.append((refSeq, sp))
    return refSeq_sp_list

#==============================================================================
# Proteins interactions
#==============================================================================

class PsicquicService:
    def __init__(self, name, restUrl):
        self.name = name
        self.restUrl = restUrl

def readURL(url):
    try:
        fileHandle = urlopen(url)
        content = fileHandle.read()
        fileHandle.close()
    except IOError:
        print('Cannot open URL ' + url)
        content = ''
    return content

def readActiveServicesFromRegistry():
    registryActiveUrl = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=xml'
    content = readURL(registryActiveUrl)
    # Create the XML reader
    root = ET.fromstring(content)
    xmlns = '{http://hupo.psi.org/psicquic/registry}'
    services = []
    for service in root.findall(xmlns + 'service'):
        name = service.find(xmlns + 'name')
        restUrl = service.find(xmlns + 'restUrl')
        service = PsicquicService(name.text, restUrl.text)
        services.append(service)
    return services

def queryPsicquic(psicquicRestUrl, query, offset, maxResults):
    query = quote(query)
    psicquicUrl = psicquicRestUrl + 'query/' + query + '?firstResult=' + str(offset) + '&maxResults=' + str(maxResults) + '&format=tab27';
    psicquicResultLines = readURL(psicquicUrl)
    return psicquicResultLines 

def parse_psicquic_result(psicquicResultLines, query_alias, representative, specie, db) -> list:
    all_ppi = []
    for line in str(psicquicResultLines).split("\n"):
        my_prot = ""
        interactant = ""
        evidence = ""
        psi = ""
        mp = ""
        cols = line.split("\t")
        # Search of infos of interest on the non uniforms psicquic output query
        if len(cols) > 1:
            for entry in cols:
                if "pubmed" in entry:
                    evidence = entry
                if "psi-mi" in entry:
                    psi = entry
            alias_a = [i for i in f"{cols[0]}|{cols[2]}|{cols[4]}".split("|") if i != "-"]
            alias_b = [i for i in f"{cols[1]}|{cols[3]}|{cols[5]}".split("|") if i != "-"]
            for alias in alias_a:
                print(alias)
                if alias.split(":")[1] in query_alias: mp = "a"
            for alias in alias_b:
                if alias.split(":")[1] in query_alias: mp = "b"
            sort_order = ["uniprotkb", "refseq", "complex", "intact"]
            alias_a = [alias.split(":")[1] for db in sort_order for alias in alias_a if alias.split(":")[0] == db]
            alias_b = [alias.split(":")[1] for db in sort_order for alias in alias_b if alias.split(":")[0] == db]
            if mp == "a": 
                my_prot = alias_a[0]
                interactant = alias_b[0].split("-")[0]
            elif mp == "b":
                my_prot = alias_b[0]
                interactant = alias_a[0].split("-")[0]
            # Stock all these info
            prot_dict = {
                "refSeq" : representative,
                "uniprot" : my_prot, 
                "func" : interactant, 
                "specie": specie,
                "db" : db, 
                "evidence" : evidence, 
                "psi-mi" : psi
                }
            all_ppi.append(prot_dict)
            print('\t' + my_prot + ' interacts with ' + interactant)
    return all_ppi
        
def query_all_psicquic_db(dict_specie_geneList) -> dict:
    """
    Use the rest psicquic api to query all interactions databases
    implementation based on
    https://github.com/PSICQUIC/psicquic-solr-ws/blob/master/src/example/python/read-all-psicquic-python3.py
    """
    directory_name = f"psicquic"
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    id_listOfInteractants = {}
    for specie, gene_list in dict_specie_geneList.items():
        # Iterate on all gene (multi refseq) of a specie
        for gene in gene_list:
            for k,v in gene.items():
                refSeq_list = v
                representative = k
            interacts_list = []
            all_id = [i for i in refSeq_list]
            for i in refSeq_list:
                all_id.extend(get_alias(i))
            # Build the complex query with alias
            conv_id = " ".join(all_id)
            query = f"alias:({conv_id}) AND species:{specie}"
            # Init services
            services = readActiveServicesFromRegistry()
            interacts_list = []
            for service in services:
                psicquic_file = pathlib.Path(f"{directory_name}/{service.name}_{specie}_{representative}.txt")
                if psicquic_file.exists():
                    with open(psicquic_file, "r") as p_file:
                        psicquic_result = ""
                        for line in p_file:
                            psicquic_result += f"{line}\n"
                else:
                    time.sleep(1)
                    psicquic_result = queryPsicquic(service.restUrl, query, 0, 200) 
                    if isinstance(psicquic_result, bytes):
                        psicquic_result = psicquic_result.decode("utf-8")
                        with open(psicquic_file, "w") as p_file:
                            p_file.write(psicquic_result)
                # Get all ppi infos from the psicquic results
                interacts_list.extend(parse_psicquic_result(psicquic_result, all_id, representative, specie, service.name))
            # Add the protein entry in the dictionnary
            id_listOfInteractants[representative] = interacts_list
    # Return all the ppi infos
    return id_listOfInteractants

def get_alias(db_id) -> list:
    all_id = [db_id]
    u = UniProt(verbose=False)
    map_refSeq_acc = u.mapping(fr="P_REFSEQ_AC", to="ID", query=db_id)
    for refSeq, id_list in map_refSeq_acc.items():
        all_id.extend(id_list)
    map_refSeq_acc = u.mapping(fr="P_REFSEQ_AC", to="ACC", query=db_id)
    for refSeq, id_list in map_refSeq_acc.items():
        all_id.extend(id_list)
    map_refSeq_acc = u.mapping(fr="P_REFSEQ_AC", to="P_ENTREZGENEID", query=db_id)
    for refSeq, id_list in map_refSeq_acc.items():
        all_id.extend(id_list)
    map_refSeq_acc = u.mapping(fr="P_REFSEQ_AC", to="UNIGENE_ID", query=db_id)
    for refSeq, id_list in map_refSeq_acc.items():
        all_id.extend(id_list)
    all_id = [i for i in all_id if i != "null"]        
    return all_id

#==============================================================================
# Subset of gene selection functions
#==============================================================================

def get_refseq_from_fasta(fasta_file) -> list:
    """
    Extract refseq of interest from a fasta file
    fasta header : >XP024837785.1_9913
    we want all the refseq with _, so XP_024837785.1
    """
    # Init list
    refseq_interest_list = []
    # Open fasta file
    with open(fasta_file, "r") as f_file:
        # Iterate on it
        for line in f_file:
            # Select header line
            if line.startswith(">"):
                # get refseq
                refseq = line.split("_")[0].replace(">","").replace("P","P_")
                # Append it to the list
                refseq_interest_list.append(refseq)
    # Return the refseq list
    return refseq_interest_list

def select_gene(dict_specie_geneList, refseq_interest_list) -> dict:
    """
    Select the gene containing our refseq of interests
    """
    set_refseq = set(refseq_interest_list)
    # Init a selction dict
    select_dict_specie_geneList = {}
    # Iterate on the full dict
    for specie, geneList in dict_specie_geneList.items():
        # For all gene dict
        for gene_dict in geneList:
            # Get the refSeq list of this gene
            for gene, refSeq_list in gene_dict.items():
                # If this gene contains 1 of my refseq of interest
                representatives = list(set_refseq.intersection(set(refSeq_list)))
                if len(representatives) > 0:
                    # Build specific dict, representative : refseq_list
                    gene_dict = {representatives[0] : refSeq for locus, refSeq in gene_dict.items()}
                    # Add it to the selection dict
                    if specie in select_dict_specie_geneList:
                        select_dict_specie_geneList[specie].append(gene_dict)
                    else:
                        select_dict_specie_geneList[specie] = [gene_dict]
    # Return the dict of the selection
    return select_dict_specie_geneList

#==============================================================================
# Write csv output file
#==============================================================================

def write_fusephylotree_csv(dict_id_IPPList, dict_specie_geneList, refseq_interest_list, filename) -> None:
    """
    Write a csv output file that could be use as fuse-phylotree input
    ex : XP_012810820.2, P59509 | P99999 etc...
    the sequence will have attributed all the function of his gene 
    (all functions of all the refseq of a given gene)
    """
    with open(filename, "w+") as csv_file:
        for specie, geneList in dict_specie_geneList.items():
            for gene_dict in geneList:
                ipp_list = []
                for gene, refSeq_list in gene_dict.items():
                    for refSeq in refSeq_list:
                        if refSeq in refseq_interest_list:
                            refSeq_interest = refSeq
                        if refSeq in dict_id_IPPList:
                            for ipp in dict_id_IPPList[refSeq]:
                                if ipp["func"] not in ipp_list:
                                    ipp_list.append(ipp["func"])
                csv_file.write(f"{refSeq_interest},{'|'.join(ipp_list)}\n")

def write_all_info_csv(dict_uniprot_IPPList, filename) -> None:
    """
    Write an output csv file with psicquic all informations
    """
    with open(filename, "w+") as csv_file:
        csv_file.write("refSeq_A,uniprot_A,uniprot_B,specie,database,evidence(s),pubmed\n")
        for prot, ipp_list in dict_uniprot_IPPList.items():
            for ipp in ipp_list:
                csv_file.write(f"{ipp['refSeq']},{ipp['uniprot']},{ipp['func']},{ipp['specie']},{ipp['db']},{ipp['psi-mi']},{ipp['evidence']}\n")

#==============================================================================
# Argparse parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("locus_protRefSeq",
                        help = "csv file, containing proteins refSeq id of a given genomic locus (ex : locus_nbProt.csv)",
                        type=str)
    parser.add_argument("fasta",
                        help = "fasta file of proteins of interest (header ; >XP024837785.1_9913), must be in the locus refseq file. If provided, only the annotations of their gene will be searched",
                        type=str)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

if __name__ == "__main__":
    # Get argparse arguments
    args = parser()
    # Read and make dict specie : listof the gene ( a gene being a dict)
    dict_specie_geneList = read_locus_prots_csv(args.locus_protRefSeq)
    # Get the refseq of these proteins of interest
    refseq_interest_list = get_refseq_from_fasta(args.fasta)
    # Make the selection of only the gene of these refseq
    dict_specie_geneList = select_gene(dict_specie_geneList, refseq_interest_list)
    # Query for all in these query list
    dict_id_IPPList = query_all_psicquic_db(dict_specie_geneList)
    # Write outputs files
    write_fusephylotree_csv(
                dict_id_IPPList,
                dict_specie_geneList, 
                refseq_interest_list,
                f"leaf_PSICQUIC_{args.fasta.split('/')[-1].split('.')[0]}.csv"
                )
    write_all_info_csv(
                dict_id_IPPList, 
                f"ipp_PSICQUIC_all_info_{args.fasta.split('/')[-1].split('.')[0]}.csv"
                )

