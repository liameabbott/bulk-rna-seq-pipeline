#!/bin/bash

# this is the "PRJ" ID, e.g. "./copy-sra-fastqs.sh PRJNA327986"
sra_id=$1
sra_dir="/efs/bioinformatics/sra/${sra_id}"
fastq_dir="${sra_dir}/fastqs"

if [ ! -d "${fastq_dir}" ]; then
  mkdir -p "${fastq_dir}";
fi

cd "${fastq_dir}"

/shared/software/bin/esearch -db sra -query $sra_id | \
/shared/software/bin/efetch -format runinfo | \
tr ',' '\t' > "${sra_id}.accessions_metadata.tsv"

grep -E "^SRR" "${sra_id}.accessions_metadata.tsv" | \
cut -f1 | \
xargs -L 1 -I {} sh -c '
    /shared/software/bin/prefetch --max-size u --progress {};
    /shared/software/bin/vdb-validate {} &> {}.prefetch_check;
    if grep -q "err" {}.prefetch_check; then
        touch {}.prefetch_failed;
    else
        touch {}.prefetch_validated;
        if ls {}*fastq.gz >/dev/null 2>&1; then
	    :
	else
	    /shared/software/bin/fasterq-dump --force --split-3 --skip-technical --print-read-nr {};
	    gzip {};
	fi
    fi
    rm {}.prefetch_check;
'

