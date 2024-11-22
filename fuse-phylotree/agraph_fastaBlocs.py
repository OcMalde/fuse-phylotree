#!/usr/bin/env python3

"""
Extract blocks from a yaml file (agraph output of paloma-2) and build a module fasta file for valide blocks.
"""


import argparse
import yaml
import os
from pathlib import Path
import random

#==============================================================================
# Extract block sequence, (from A Belcour, slightly modified)
#==============================================================================

def extract_block_sequence(input_file):
    """
    Extract for each PLMA block its sequence found in the input sequences.
    
    Parameters
    ----------
    input_file : str
        Name of a file in yaml format
        
    Returns
    -------
    dict_bloc_infos : dict
        Dictionary with bloc as key, bloc informations as value
    seqNumber_seqInfos : dict
        Dictionary with sequence number as key, sequence informations as value
    """
    # Read the yaml file
    with open(input_file) as i_file:
        agraph_data = yaml.load(i_file, Loader=yaml.FullLoader)
    # Get sequences infos (number / name / sequence)
    seqNumber_seqInfos = {}
    for seq_number, dict_infos in agraph_data["Sequences"].items():
        seq_name = dict_infos["id"]
        seq = dict_infos["seq"]
        seqNumber_seqInfos[seq_number] = (seq_name, seq)
    # Get all blocs infos / and get their sequences from their positions
    dict_bloc_infos = {}
    for zone, dict_zone in agraph_data["Aligned sites"].items():
        for bloc, list_bloc in dict_zone.items():
            blocs_infos = {}
            if len(list_bloc) >= 5 :
                for res in list_bloc:
                    for seq_number, res_pos in res.items():
                        res_pos = res_pos[0]
                        if seq_number not in blocs_infos:
                            blocs_infos[seq_number] = {"start": int(res_pos), "end": int(res_pos), "seq" : ""}
                        else:
                            if int(res_pos) > blocs_infos[seq_number]["end"]:
                                blocs_infos[seq_number]["end"] = res_pos
                            elif int(res_pos) < blocs_infos[seq_number]["start"]:
                                blocs_infos[seq_number]["start"] = res_pos
                for seq_number in blocs_infos:
                    blocs_infos[seq_number]["seq"] = seqNumber_seqInfos[seq_number][1][blocs_infos[seq_number]["start"]-1:blocs_infos[seq_number]["end"]]
                # TODO ; keep only if len >= 2
                if len(blocs_infos[seq_number]["seq"]) > 1:
                    dict_bloc_infos[f"Z{zone}B{bloc}"] = blocs_infos
    return dict_bloc_infos, seqNumber_seqInfos

#==============================================================================
# Build bloc fasta file
#==============================================================================

def writeFastas(dict_bloc_infos, seqNumber_seqInfos, fileName, directory) -> None:
    """
    Build a fasta file for each paloma bloc
    
    Parameters
    ----------
    dict_bloc_infos : dict
        Dictionary with bloc as key, bloc informations as value
    seqNumber_seqInfos : dict
        Dictionary with sequence number as key, sequence informations as value
    fileName : str
        Name of a file in yaml format
    directory : str
        Name of the directory where bloc/modules fasta files will be written    
    """
    family = Path(fileName).stem.split("_")[0]
    if not os.path.exists(directory):
        os.makedirs(directory)
    for block, infos in dict_bloc_infos.items():
        with open(f"{directory}/{block}.fasta","w") as f_file:
            for seq_id, seq_info in infos.items():
                seq_name = seqNumber_seqInfos[seq_id][0]
                seq = seq_info["seq"]
                start = seq_info["start"]
                end = seq_info["end"]
                f_file.write(f">{block}|{start}|{end}|_{seq_name.split(' ')[0].split('_')[0]}_{family}\n{seq}\n")
        # If there is only 2 sequences in the blocs, make a trivial tree file
        if len(infos) == 2:
            names = []
            for seq_id, seq_info in infos.items():
                seq_name = seqNumber_seqInfos[seq_id][0]
                start = seq_info["start"]
                end = seq_info["end"]
                names.append(f"{block}|{start}|{end}|_{seq_name.split(' ')[0].split('_')[0]}_{family}")
            with open(f"{directory}/{block}.tree","w") as t_file:
                t_file.write(f"({names[0]},{names[1]});")

#==============================================================================
# Treat a dot file and write the modules files on a specific given directory
#==============================================================================

def make_module_directory(in_file, module_directory) -> None:
    """
    Extract blocs from a yaml file and write the corresponding fasta files
    
    Parameters
    ----------
    in_file : str
        Name of a file in yaml format
    module_directory : str
        Name of the directory where bloc/modules fasta files will be written  
    """
    blocks_list, num_name = extract_block_sequence(in_file)
    writeFastas(blocks_list, num_name, in_file, module_directory)

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("agraph_file", help = "yaml (.agraph) result file")
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    blocks_list, num_name = extract_block_sequence(args.agraph_file)
    writeFastas(blocks_list, num_name, args.agraph_file, f"blocs_fasta")



if __name__ == '__main__':
	main()
