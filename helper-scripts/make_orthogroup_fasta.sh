#!/bin/bash
#SBATCH --job-name=advanced
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --output=sortie_output.out

# Enable strict error handling
set -euo pipefail

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <id_seq_file> <name> <docker_image|singularity_image_path>"
    exit 1
fi

# Convert input arguments to absolute paths
id_seq_file=$(realpath "$1")        # File with IDs of proteins
name="$2"                           # Name of the protein family (e.g., "fibulin")
container="$3"                      # Docker image name or path to Singularity image

# Determine container type (Docker or Singularity)
if [[ "$container" == *.sif ]] || [[ "$container" == /* ]]; then
    container_type="singularity"
    container_path=$(realpath "$container")
else
    container_type="docker"
    container_path="$container"
fi

# Current working directory
workplace=$(pwd)

####################################
####### Make FASTA file ###########
####################################
mkdir -p run_orthogroup
cd run_orthogroup

# Commands for FASTA creation
cmd1="python3 /fuse-phylotree/myOrthogroups_fasta.py --download /data_9sp/Orthogroups.tsv ${id_seq_file} /data_9sp/assocF_taxid_spName.csv"
cmd2="python3 /fuse-phylotree/gff_regroup_iso_locus.py --fasta_directory orthogroups_*/ --assoc_file /data_9sp/assocF_taxid_dbnt.csv --gff_directory /data_9sp/gff/"

# Execute commands based on container type
if [ "$container_type" == "docker" ]; then
    docker run --rm -v "${workplace}:${workplace}" -w "$(pwd)" ${container_path} bash -c "${cmd1}"
    docker run --rm -v "${workplace}:${workplace}" -w "$(pwd)" ${container_path} bash -c "${cmd2}"
else
    singularity exec "${container_path}" ${cmd1}
    singularity exec "${container_path}" ${cmd2}
fi

# Optional: Create taxid list
grep ">" ${workplace}/run_orthogroup/orthogroups_*/isoforms_per_locus/longest_isoform.fasta \
    | sed -E 's/^>[^_]+_//' | sort -u > ${workplace}/taxid_list.txt

cd "$workplace"

####################################
########### Copy FASTA #############
####################################
# Get generated FASTA file
path_fasta=$(find "${workplace}/run_orthogroup" -name longest_isoform.fasta | head -n 1)
num_seqs=$(grep -c ">" "$path_fasta")
output_fasta="${name}_${num_seqs}.fasta"

cp "$path_fasta" "${workplace}/${output_fasta}"

echo "FASTA file created: ${output_fasta}"

