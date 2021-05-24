#!/bin/bash

# memory environment variable for Nextflow
export NXF_OPTS="-Xms1g -Xmx4g"

# script takes PRJ ID as an argument
sra_id=$1
fastq_dir="/efs/bioinformatics/sra/${sra_id}/fastqs"
output_dir="/efs/bioinformatics/sra/${sra_id}/processed-data"

if [ ! -d "${output_dir}" ]; then
  mkdir $output_dir;
fi

cd /shared/bioinformatics/bulk-rna-seq-pipeline

/shared/software/nextflow run main.nf \
  -latest \
  -profile conda,slurm \
  --dataset_name "${sra_id}" \
  --fastq_directory "${fastq_dir}" \
  --output_directory "${output_dir}" 

