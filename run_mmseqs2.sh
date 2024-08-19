#!/bin/bash

#SBATCH --job-name=mmseqs2
#
#SBATCH --time=3-00:00:00 # days-hh:mm:ss
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=41
#SBATCH --mem-per-cpu=2048 # 4GB
#
#SBATCH --mail-user=mikhail.v.fofanov@gmail.com
#SBATCH --mail-type=ALL

#
#SBATCH --output=result_mmseqs2.%j.txt

date;hostname;pwd

#. $CONDA_ROOT/etc/profile.d/conda.sh
source /home/zo49sog/miniconda3/etc/profile.d/conda.sh

conda activate

#module load tools/python/3.8

cd /home/zo49sog/crassvirales/phylomes/

threads=40

prots_faa="/home/zo49sog/crassvirales/phylomes/crassvirales_refseq/crassvirales.faa"

analysis_id='crassvirales_refseq'

tsv="results/${analysis_id}/6_protein_clustering/table_clustering.tsv"
tmp_dir="results/${analysis_id}/6_protein_clustering/tmp"
prots_db_dir="results/${analysis_id}/6_protein_clustering/db"
prots_db="results/${analysis_id}/6_protein_clustering/db/all_proteins"
out_prefix="results/${analysis_id}/6_protein_clustering/clustering"

logs_dir="logs/${analysis_id}/protein_clustering/"
log_db="${logs_dir}/db.log"
log_clust="${logs_dir}/clustering.log"
log_tsv="${logs_dir}/to_table.log"

mkdir -p ${prots_db_dir}
mkdir -p ${logs_dir}

mmseqs createdb ${prots_faa} ${prots_db} >> ${log_db}
mmseqs cluster ${prots_db} ${out_prefix} ${tmp_dir}\
	-c 0.8 --min-seq-id 0 -s 7.5 -e '1.000E-03' --cluster-steps 4 --threads ${threads}\
	--cluster-reassign >> ${log_clust}
mmseqs createtsv ${prots_db} ${prots_db} ${out_prefix}\
	${tsv} >> ${log_tsv}

conda deactivate

date
