#!/bin/bash

#SBATCH --job-name=advanced
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --output=sortie_output.out

# Activer la gestion des erreurs
set -euo pipefail

# Path to singularity image 
path_image="fuse_phylotree.sif"
workplace=$(pwd)

####### File path 

# File with id of protein of interest
id_seq_file=${workplace}/human_fibulins_ids.txt

# File with id of the outgroup
id_outgroup_file=${workplace}/human_outgroup_ids.txt

####### Variable 

# Name of protein of interest 
name="fibulin"

##################################################
############### Make fasta file ##################
############# of protein interest ################
##################################################

# Make a folder
mkdir run_orthogroup
cd run_orthogroup

# Command to run
cmd1="python3 /fuse-phylotree/myOrthogroups_fasta.py --download /data_9sp/Orthogroups.tsv ${id_seq_file} /data_9sp/assocF_taxid_spName.csv"

echo $cmd1
cmd2="python3 /fuse-phylotree/gff_regroup_iso_locus.py --fasta_directory orthogroups_*/ --assoc_file /data_9sp/assocF_taxid_dbnt.csv --gff_directory /data_9sp/gff/"

# Run the script
singularity exec "${path_image}" ${cmd1}

singularity exec "${path_image}" ${cmd2}

# Make taxid list from file fasta 

grep ">" ${workplace}/run_orthogroup/orthogroups_*/isoforms_per_locus/longest_isoform.fasta | sed -E 's/^>[^_]+_//' | sort -u > ${workplace}/taxid_list.txt
taxid_list="${workplace}/taxid_list.txt"

cd "$workplace"

##################################################
############### Make fasta file ##################
################# of outgroup ####################
##################################################
mkdir outgroup_fasta
cd outgroup_fasta

# Command to run
cmd1="python3 /fuse-phylotree/myOrthogroups_fasta.py --download /data_9sp/Orthogroups.tsv "${id_outgroup_file}" /data_9sp/assocF_taxid_spName.csv"

cmd2="python3 /fuse-phylotree/gff_regroup_iso_locus.py --fasta_directory orthogroups_*/ --assoc_file /data_9sp/assocF_taxid_dbnt.csv --gff_directory /data_9sp/gff/"

# Run the script
singularity exec "${path_image}" ${cmd1}

singularity exec "${path_image}" ${cmd2}

cd "$workplace"

##################################################
###### Make species tree from a taxid list #######
##################################################

# Make a folder 

mkdir dir_tree
cd dir_tree

# Command to run 
cmd="python3 /fuse-phylotree/species_phylo.py ${taxid_list}"

# Run the script 
singularity exec "${path_image}" ${cmd} > out_species_tree.txt

##################################################
############## Make gene tree ####################
###### from a fasta file and a species tree ######
##################################################

# Path to the input file
path_file_fasta="${workplace}/run_orthogroup/orthogroups_*/isoforms_per_locus/longest_isoform.fasta"

# Copy input file in worplace
number_sequence=$(grep -c ">" "$path_file_fasta")
name_fasta_file="${name}_${nummber_sequence}.fasta"
cp $path_file_fasta ${workplace}/${name_fasta_file}

# Path to outgroup fasta
path_output_fasta="${workplace}/outgroup_fasta/orthogroups_*/isoforms_per_locus/longest_isoform.fasta"

# Make fasta file with interest protein and outgroup
file_out_fbln="${name}_outgroup.fasta"
cp ${workplace}/${name_fasta_file} $file_out_fbln
cat $path_output_fasta >> $file_out_fbln

species="*_species.tree"

# Command to run
cmd="python3 /fuse-phylotree/gene_phylo.py $file_out_fbln $species"

# Run the scirpt 
singularity exec "${path_image}" ${cmd}

cp "gene_phylo_dir_*/outgroup.tree" "${workplace}/no_root_tree.tree"

cd "$workplace"
