#!/bin/bash
#SBATCH --job-name=advanced
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --output=sortie_output.out

# Enable error handling
set -euo pipefail

# Check for correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <id_seq_file> <id_outgroup_file> <name> <docker_image>"
    exit 1
fi

# Convert input file paths to absolute paths
id_seq_file=$(realpath "$1")          # File with IDs of proteins of interest (in-group)
id_outgroup_file=$(realpath "$2")       # File with IDs of outgroup proteins
name="$3"                             # Name of the protein family (e.g., "fibulin")
docker_image="$4"                     # Docker image name (e.g., ghcr.io/ocmalde/fuse-phylotree:test)

# Get the current working directory
workplace=$(pwd)

##################################################
############### Make FASTA file ##################
######### for protein of interest ##############
##################################################
mkdir run_orthogroup
cd run_orthogroup

# Command to run for in-group proteins
cmd1="python3 /fuse-phylotree/myOrthogroups_fasta.py --download /data_9sp/Orthogroups.tsv ${id_seq_file} /data_9sp/assocF_taxid_spName.csv"
echo $cmd1
cmd2="python3 /fuse-phylotree/gff_regroup_iso_locus.py --fasta_directory orthogroups_*/ --assoc_file /data_9sp/assocF_taxid_dbnt.csv --gff_directory /data_9sp/gff/"

docker run --rm -v "${workplace}:${workplace}" -w "$(pwd)" ${docker_image} bash -c "${cmd1}"
docker run --rm -v "${workplace}:${workplace}" -w "$(pwd)" ${docker_image} bash -c "${cmd2}"

# (Optional) Create a taxid list if needed later
grep ">" ${workplace}/run_orthogroup/orthogroups_*/isoforms_per_locus/longest_isoform.fasta | sed -E 's/^>[^_]+_//' | sort -u > ${workplace}/taxid_list.txt

cd "$workplace"

##################################################
############### Make FASTA file ##################
################# for outgroup ###################
##################################################
mkdir outgroup_fasta
cd outgroup_fasta

# Command to run for outgroup proteins
cmd1="python3 /fuse-phylotree/myOrthogroups_fasta.py --download /data_9sp/Orthogroups.tsv ${id_outgroup_file} /data_9sp/assocF_taxid_spName.csv"
cmd2="python3 /fuse-phylotree/gff_regroup_iso_locus.py --fasta_directory orthogroups_*/ --assoc_file /data_9sp/assocF_taxid_dbnt.csv --gff_directory /data_9sp/gff/"

docker run --rm -v "${workplace}:${workplace}" -w "$(pwd)" ${docker_image} bash -c "${cmd1}"
docker run --rm -v "${workplace}:${workplace}" -w "$(pwd)" ${docker_image} bash -c "${cmd2}"

cd "$workplace"

##################################################
########### Create Combined FASTA File ###########
##################################################
# Get the in-group FASTA file
path_in_fasta=$(find "${workplace}/run_orthogroup" -name longest_isoform.fasta | head -n 1)
num_in=$(grep -c ">" "$path_in_fasta")
in_fasta="${name}_${num_in}.fasta"
cp "$path_in_fasta" "${workplace}/${in_fasta}"

# Get the outgroup FASTA file
path_out_fasta=$(find "${workplace}/outgroup_fasta" -name longest_isoform.fasta | head -n 1)
num_out=$(grep -c ">" "$path_out_fasta")
out_fasta="outgroup_${name}_${num_out}.fasta"
cp "$path_out_fasta" "${workplace}/${out_fasta}"

# Create the combined FASTA file by concatenating the two files
combined_num=$(( num_in + num_out ))
combined_fasta="${name}_${combined_num}.fasta"
cp "${workplace}/${in_fasta}" "${workplace}/${combined_fasta}"
cat "${workplace}/${out_fasta}" >> "${workplace}/${combined_fasta}"

echo "Combined FASTA file created: ${combined_fasta} | In-group: ${in_fasta} | Out-group: ${out_fasta}"
