#!/bin/bash

#SBATCH --job-name=snakemake_tree_analysis
#SBATCH --output=/home/zo49sog/crassvirales/veo_phage_genome_analysis/slurm_logs/result_%x.%j.txt
#SBATCH --time=3:00:00  # This will be overwritten by the script
#SBATCH --partition=short  # This will be overwritten by the script
#SBATCH --nodes=1  # This will be overwritten by the script
#SBATCH --ntasks=1  # This will be overwritten by the script
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB

date; hostname; pwd

# Load Snakemake environment or module
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh

working_dir="/mnt/c/crassvirales/Bas_phages_large/Bas_phages/5_nr_screening/4_merged_ncbi_crassvirales/2_trees_leaves/tree_analysis_test"

snakefile="/mnt/c/crassvirales/phylomes/scripts/tree_analysis/Snakefile"
jobs=5000
latency_wait=20
slurm_output_dir="${working_dir}/results"
slurm_time="3:00:00"
slurm_partition="short"
slurm_nodes=1
slurm_ntasks=1

# Change to the working directory
cd ${working_dir}

snakemake --snakefile "${snakefile}" --jobs ${jobs} --cluster "sbatch --output=${slurm_output_dir}/{rule}_%x_%j.out --time=${slurm_time} --partition=${slurm_partition} --nodes=${slurm_nodes} --ntasks=${slurm_ntasks} --cpus-per-task={threads}" --latency-wait ${latency_wait}

# Run Snakemake with SLURM cluster submission
#snakemake --snakefile "${snakefile}" --jobs ${jobs} --cluster "sbatch --output=${slurm_output_dir}/{rule}_%x_%j.out --time=${slurm_time} --partition=${slurm_partition} --nodes=${slurm_nodes} --ntasks=${slurm_ntasks} --cpus-per-task={threads} --mem={resources.mem_mb}MB" --latency-wait ${latency_wait}

date