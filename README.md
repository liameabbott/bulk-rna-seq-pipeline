# bulk-rna-seq-pipeline

This Nextflow pipeline takes raw FASTQ reads from bulk RNA-seq experiments and produces transcript- and gene-by-sample matrices of read counts to be used in downstream analysis activities, such as differential gene expression studies.

### Software Requirements

* [Nextflow](https://www.nextflow.io/)
* [SRA Toolkit](https://github.com/ncbi/sra-tools)

All other software requirements are managed by the conda environment (`conda/environment.yml`) used by the pipeline.

### Example Usage

Follow these steps to download a dataset from SRA, process it with bulk RNA-seq pipeline, and copy the processed files to S3. 

1. SSH in to the Slurm cluster head node:
```
ssh -i <path-to-your-pem-key> <username>@mri-hpc-slurm.gatesmri.site
```

2. Navigate to the cloned copy of this repository on the `/shared` drive:
```
cd /shared/bioinformatics/bulk-rna-seq-pipeline
```

3. Start a `screen` or `tmux` session to decouple pipeline processes from the terminal shell process:
```
screen -a
```

4. (optional) If not already downloaded, use the script `scripts/copy-sra-fastqs.sh` to download SRA files for a given project (e.g. `PRJNA327986`) from the SRA repository and convert them to FASTQs:
```
./scripts/copy-sra-fastqs.sh PRJNA327986
```
   This script will download SRA files to `/efs/bioinformatics/sra/cache` and converted FASTQ files will be written to `/efs/bioinformatics/sra/<project-id>/fastqs/` (e.g. `/efs/bioinformatics/sra/PRJNA327986/fastqs/`).
   
   Note that this step does not need to be run from the cluster head node -- it could also be run from the `bio-dev` server.

5. Make sure you are logged in to the Slurm cluster head node and that you are in a `screen` or `tmux` session. Then, run the pipeline with default reference database settings (Ensembl, release 102, Homo sapiens, GRCh38 reference files) for a given project using the script `scripts/run-bulk-rna-seq-pipeline.sh`:
```
./scripts/run-bulk-rna-seq-pipeline.sh PRJNA327986
```

   The full Nextflow command run by this runner script is:
   ```
   /shared/software/nextflow run main.nf -latest -profile conda,slurm \
   --dataset_name PRJNA327986 \
   --fastq_directory /efs/bioinformatics/sra/PRJNA327986/fastqs \
   --output_directory /efs/bioinformatics/sra/PRJNA327986/processed-data
   ```

   The `-profile conda,slurm` argument tells the pipeline to use the conda environment defined by `conda/environment.yml` and to distribute the work over cores in the Slurm cluster.
    
