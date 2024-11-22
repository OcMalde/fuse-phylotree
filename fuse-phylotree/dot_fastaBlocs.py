#!/usr/bin/env python3

"""
@author: odennler
"""


import argparse
import functools
import os
import random
import shutil
from pathlib import Path

#==============================================================================
# Extract block sequence, (from A Belcour, slightly modified)
#==============================================================================

def extract_block_sequence(plma_input_file, blockList_sequence_num_ids, seqNumber_seqName):
    """
    Extract for each PLMA block its sequence found in the input sequences.
    """
    concatDotFile = functools.reduce(lambda line1, line2: line1.strip()+line2.strip(), open(plma_input_file,'r'))
    # Create a list with each subgraph in it by splitting on subgraph.
    # Extract the first cluster (containing sequence numbers and their name)
    subgraph_seqName = concatDotFile.split('subgraph cluster_')[1]
    subgraph_clust = subgraph_seqName.split(';"')[1:]
    for clust in subgraph_clust:
        clust_data = clust.split('[label = "')[1].split('"]')[0]
        seqNum = clust_data.split(": ")[0]
        seqName = clust_data.split(": ")[1]
        seqNumber_seqName[seqNum] = seqName
    # Extract from the third cluster (the first cluster containing PLMA blocks.)
    subgraph_clusters = concatDotFile.split('subgraph cluster_')[3:]
    # For each subgraph of the dot file (=PLMA block after block 2).
    # Extract the sequence_number and the sequence_identifier.
    for subCluster in subgraph_clusters:
        # Remove the data at the end of the dot file.
        subgraph_data = subCluster.split('}')[0]
        # Extract the number of the block which just after 'subgraph cluster_'.
        #block_number = subgraph_data.split('{')[0]
        # Then creates a list with all the line inside the {} of the 'subgraph cluster_'.
        subgraph_lines = subgraph_data.split(';"')[1:]
        seq_num_ids = {}
        for subgraph_line in subgraph_lines:
            # Extract the sequence number (first number in the "(SEQNUM,START,END)" [label = "X: SEQID"];)
            num_sequence = subgraph_line.split('" [label = "')[0].strip('"').split(',')[0].strip('(')
            # Extract the block start position.
            start_block_in_seq = subgraph_line.split('" [label = "')[0].strip('"').split(',')[1].strip('(')
            # Extract the block end position (by extract the size and add it to the start)
            end_block_in_seq = subgraph_line.split('" [label = "')[0].strip('"').split(',')[2].strip(')')
            end_block_in_seq = int(start_block_in_seq) + int(end_block_in_seq) - 1
            # Extract the block sequence.
            block_sequence = subgraph_line.split('" [label = "')[1].split('"]')[0]
            seq_num_ids[num_sequence] = (block_sequence, start_block_in_seq, end_block_in_seq)
        blockList_sequence_num_ids.append(seq_num_ids)
    return blockList_sequence_num_ids, seqNumber_seqName

#==============================================================================
# Build bloc fasta file
#==============================================================================

def writeFastas(blockList_sequence_num_ids, proteinNameDict, thres, fileName, directory) -> None:
    """
    Build a fasta file for each paloma bloc
    """
    family = Path(fileName).stem.split("_")[0]
    if not os.path.exists(directory):
        os.makedirs(directory)
    sequence = ""
    block = 0
    for seq_num_ids in blockList_sequence_num_ids:
        block += 1
        # If the sequence of the bloc is >= threshold, write the blocs in a multi fasta file
        if len(random.choice(list(seq_num_ids.values()))[0]) >= int(thres):
            with open(f"{directory}/B{block}.fasta","w") as f_file:
                bloc_info = []
                for seq_ids, seq_info in seq_num_ids.items():
                    sequence = seq_info[0]
                    start = int(seq_info[1])
                    end = int(seq_info[2])
                    seq_name = proteinNameDict[seq_ids]
                    f_file.write(f">B{block}|{start}|{end}|_{seq_name.split(' ')[0].split('_')[0]}_{family}\n{sequence}\n")
            # If there is only 2 sequences in the blocs, make a trivial tree file
            if len(seq_num_ids) == 2:
                names = []
                for seq_ids, seq_info in seq_num_ids.items():
                    sequence = seq_info[0]
                    start = int(seq_info[1])
                    end = int(seq_info[2])
                    seq_name = proteinNameDict[seq_ids]
                    names.append(f"B{block}|{start}|{end}|_{seq_name.split(' ')[0].split('_')[0]}_{family}")
                with open(f"{directory}/B{block}.tree","w") as t_file:
                    t_file.write(f"({names[0]},{names[1]});")

#==============================================================================
# Treat a dot file and write the modules files on a specific given directory
#==============================================================================

def make_module_directory(dot_file, thres, module_directory) -> None:
    blocks_list, num_name = [], {}
    blocks_list, num_name = extract_block_sequence(dot_file, blocks_list, num_name)
    writeFastas(blocks_list, num_name, thres, dot_file, module_directory)

#==============================================================================
# Parser
#==============================================================================

def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("dot_file", help = "PLMA DOT (.dot) result file")
    parser.add_argument("thres", help = "Threshold; minimum len of a bloc")
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    args = parser()
    blocks_list, num_name = [], {}
    blocks_list, num_name = extract_block_sequence(args.dot_file, blocks_list, num_name)
    writeFastas(blocks_list, num_name, args.thres, args.dot_file, f"blocs_fasta_t{thres}")



if __name__ == '__main__':
	main()
