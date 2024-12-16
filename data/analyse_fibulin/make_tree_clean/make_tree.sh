#!/bin/bash

#SBATCH --job-name=tree
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --output=sortie_output.out

file_fasta="fibuline_outgroup_clean.fasta"
species="8_species.tree"
image_path="/home/genouest/inra_umr1349/echenel/to_github/fuse_phylotree.sif"

cmd="python3 /fuse-phylotree/gene_phylo.py $file_fasta $species"

singularity exec ${image_path} ${cmd}

