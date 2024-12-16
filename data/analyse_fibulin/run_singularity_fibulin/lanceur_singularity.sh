#!/bin/bash

#SBATCH --job-name=lanceur
#SBATCH --cpus-per-task=20
#SBATCH --mem=50G
#SBATCH --output=sortie_lanceur.out

date
image_path=/home/genouest/inra_umr1349/echenel/fuse_phylotree.sif

file_fasta="fibulin_59.fasta"

file_ppi="ppi_shared.csv"
name=$(basename "$file_fasta" | cut -d. -f1)

gene_tree="fibulin_tree_root.tree"

cmd2="python3 /fuse-phylotree/fuse-phylotree.py --output_directory dir_fibuline_phylocharmod --gene_tree $gene_tree $file_fasta $file_ppi"

sortir="sortie_fibulin_analyse.txt"

echo "Path to output : "
echo $output_dir

echo "Begin with the singularity image"
singularity exec ${image_path} ${cmd2} > $sortir

echo "Analyse finished" 
date


