#!/bin/bash

#SBATCH --job-name=lanceur
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --output=sortie_lanceur.out

# File paths
image_path="fuse-phylotree.sif"
file_fasta="fibulin_59.fasta"
file_annotation="ppi_shared.csv"
gene_tree="fibulin_tree_root.tree"
output_dir="dir_analysis_fuse_phylotree"

# Command to launch tool
cmd="python3 /fuse-phylotree/fuse-phylotree.py --output_directory ${output_dir} $file_fasta $file_annotation $gene_tree"

# Start analysis
singularity exec ${image_path} ${cmd}


