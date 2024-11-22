#!/bin/python3

import argparse
import random
import os

#==============================================================================
# Sequence class
#==============================================================================

class sequence(object):
    
    def __init__(self, header, seq):
        # In case of ncbi id, replace XP_ , NP_ by XP, NP
        header = header.replace("P_","P").replace("\n","")
        self.header = header
        self.header.split("_")[0]
        self.seq = seq.replace("\n","")

    def write(self) -> str:
        to_write = (
                str(self.header)
                + "\n"
                + str(self.seq)
                + "\n"
                )
        return to_write

#==============================================================================
# Parse multiple fasta file
#==============================================================================

def parseFasta(fileName) -> list:
    sequence_list = []
    with open(fileName, "r") as fastaFile:
        current_sequence = ""
        for line in fastaFile:
            if line.startswith(">"):
                if current_sequence != "":
                    sequence_list.append(
                            sequence(header = header, seq = current_sequence)
                            )
                    current_sequence = ""
                header = line
            else:
                current_sequence += line
        sequence_list.append(
                        sequence(header = header, seq = current_sequence)
                            )
    return sequence_list

#==============================================================================
# Reechantillonage by addition from another pull of sequence
#==============================================================================

def addR_from_pull(allFasta_list, seqPull, n_seq, n_iter) -> list:
    """
    Add random n_seq from seqPull, n_iter (= n_times)
    (= n_iter differents reechantillonage)
    And return a list with the n_iter corresponding protein list
    """
    proteinList_list = []
    # For the n_iterations, add n_seq
    for n_i in range(n_iter):
        # For the actual reech, init protein list with the initial prot list
        reech_proteinList = [i for i in allFasta_list]
        # Take n_seq random sequences from seqPull
        for n_s in range(n_seq):
            # Pick a random seq from seqPull
            # That is not already in our sequence set
            not_fund = False
            while not_fund == False:
                pick = random.choice(seqPull)
                if pick.seq not in [s.seq for s in reech_proteinList]:
                    not_fund = True
            # If not in, add the random to the actual reech
            reech_proteinList.append(pick)
        # Now that this reechentillonage is completed, add it to the list of list
        proteinList_list.append(reech_proteinList)
    return proteinList_list

#==============================================================================
# Write the all the partials fasta into files (in out directory)
#==============================================================================

def writte_all(list_Partiallist, output_directory) -> None:
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    index = 1
    for partialList in list_Partiallist:
        fileName = (str(output_directory).replace("/","")
                + "/"
                + str(index)
                + ".fasta")
        write_fasta(partialList, fileName)
        index += 1

def write_fasta(seq_list, fileName) -> None:
    with open(fileName, "w+") as fastaFile:
        [fastaFile.write(sequences.write()) 
                for sequences in seq_list]

#==============================================================================
# Parser
#==============================================================================

def parser():
    """
    Argparse parser function
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("multi_fasta",
            help = "Multiple fasta file, containing all the sequences of interest",
                        type=str)
    parser.add_argument("adit_pull",
            help = "Multiple fasta file, containing a pull of sequences to add",
                        type=str)
    parser.add_argument("n_seq",
            help = "Number of sequence to add in each reechantillonnage",
                        type=int)
    parser.add_argument("n_iter",
            help = "Number of reechantillonnage to do",
                        type=int)
    return parser.parse_args()

#==============================================================================
# Main
#==============================================================================

def main():
    """
    Main function
    """
    # Parserser call
    args = parser()
    n_seq, n_iter = args.n_seq, args.n_iter
    # Read an parse fasta into a proteins object list 
    allFasta_list = parseFasta(args.multi_fasta)
    # Read an parse fasta of aditional sequence pull
    seqPull = parseFasta(args.adit_pull)
    # Add nb random sequences of aditional pull, n times
    proteinList_list = addR_from_pull(allFasta_list, seqPull, n_seq, n_iter)
    # Write all the fastas (1 fasta per protein list in the proteinList_list)
    writte_all(proteinList_list, f"reech_{args.multi_fasta}_S{n_seq}_I{n_iter}")
    



if __name__ == '__main__':
    main()


