# bulk-rna-seq-pipeline

This Nextflow pipeline takes raw FASTQ reads from bulk RNA-seq experiments and produces transcript- and gene-by-sample matrices of read counts to be used in downstream analysis activities, such as differential gene expression studies.

### Software Requirements

* [Nextflow](https://www.nextflow.io/)
* [SRA Toolkit](https://github.com/ncbi/sra-tools)

All other software requirements are managed by the conda environment (`conda/environment.yml`) used by the pipeline.

### Downloading a dataset from SRA

Use the script `scripts/copy-sra-fastqs.sh` to download raw read data from SRA and convert it to FASTQs. These FASTQs are the input to the pipeline.

The script takes an SRA project's `PRJ` ID as an argument and downloads all samples associated with that project to `/efs/bioinformatics/sra/<PRJ-ID>/`.

#### Example Usage

`./scripts/copy-sra-fastqs.sh PRJNA327986`
