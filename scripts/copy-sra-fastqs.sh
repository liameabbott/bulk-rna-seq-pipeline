#!/bin/bash

# this is the "PRJ" ID, e.g. "./copy-sra-fastqs.sh PRJNA327986"
sra_id=$1
sra_dir="/efs/bioinformatics/sra/${sra_id}"

if [ ! -d "${sra_dir}" ]; then
  mkdir -p "${sra_dir}";
fi

cd "${sra_dir}"

/shared/software/bin/esearch -db sra -query $sra_id | \
/shared/software/bin/efetch -format runinfo | \
tr ',' '\t' > "${sra_id}.accessions_metadata.tsv"

grep -E "^SRR" "${sra_id}.accessions_metadata.tsv" | \
cut -f1 | \
xargs -L 1 -I {} sh -c '
    prefetch --max-size u --progress {};
    vdb-validate {} &> {}.prefetch_validation;
    if grep -q "err" {}.prefetch_validation; then
        touch {}.prefetch_failed;
    else
        touch {}.prefetch_validated;
        fasterq-dump --force --split-3 --skip-technical --print-read-nr {};
    fi
    rm {}.prefetch_validation;
' && \
parallel gzip ::: *.fastq

