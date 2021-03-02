
process clean_single_reads_fastp {
    label "lrg"
    publishDir "${params.output_directory}/${id}/qc", \
        mode: "copy", overwrite: true, \
        pattern: "*.fastp.{json,html}"

    input:
    tuple val(id), path(read)

    output:
    tuple val(id), path("*fastp*")

    """
    fastp \
    --in1 ${read} \
    --out1 ${id}.fastp.fastq.gz \
    --json ${id}.fastp.json \
    --html ${id}.fastp.html \
    --thread ${task.cpus}   
    """
}

process clean_paired_reads_fastp {
    label "lrg"
    publishDir "${params.output_directory}/${id}/qc", \
        mode: "copy", overwrite: true, \
        pattern: "*.fastp.{json,html}"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("*fastp*")

    """
    fastp \
    --in1 ${reads[0]} \
    --out1 ${id}_1.fastp.fastq.gz \
    --in2 ${reads[1]} \
    --out2 ${id}_2.fastp.fastq.gz \
    --detect_adapter_for_pe \
    --json ${id}.fastp.json \
    --html ${id}.fastp.html \
    --thread ${task.cpus}
    """
}

process quantify_single_reads_salmon {
    label "lrg"
    publishDir "${params.output_directory}/${id}", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "salmon" }

    input:
    tuple val(id), path(fastp_fastq)
    path(salmon_index)
    path(reference_gtf)

    output:
    tuple val(id), env(strand), path("${id}-salmon", type: "dir")

    """
    salmon quant \
    --index ${salmon_index} \
    --libType A \
    --output ${id}-salmon \
    --unmatedReads ${fastp_fastq} \
    --gcBias \
    --seqBias \
    --recoverOrphans \
    --threads ${task.cpus}

    strand=\$(
        cat ${id}-salmon/lib_format_counts.json | \
        jq '.expected_format' | \
        sed 's/"//g'
    )
    
    mv ${id}-salmon/quant.sf ${id}-salmon/${id}.quant.sf
    """
}

process quantify_paired_reads_salmon {
    label "lrg"
    publishDir "${params.output_directory}/${id}", \
        mode: "copy", overwrite: true, \
        saveAs: { filename -> "salmon" }

    input:
    tuple val(id), path(fastp_fastqs)
    path(salmon_index)
    path(reference_gtf)

    output:
    tuple val(id), env(strand), path("${id}-salmon", type: "dir")

    """
    salmon quant \
    --index ${salmon_index} \
    --libType A \
    --output ${id}-salmon \
    --mates1 ${fastp_fastqs[0]} \
    --mates2 ${fastp_fastqs[1]} \
    --gcBias \
    --seqBias \
    --recoverOrphans \
    --threads ${task.cpus}

    strand=\$(
        cat ${id}-salmon/lib_format_counts.json | \
        jq '.expected_format' | \
        sed 's/"//g'
    )
    
    mv ${id}-salmon/quant.sf ${id}-salmon/${id}.quant.sf
    """
}

process align_single_reads_star {
    label "lrg"
    publishDir "${params.output_directory}/${id}", \
        mode: "copy", overwrite: true

    input:
    tuple val(id), path(fastp_fastq)
    path(star_index)

    output:
    tuple val(id), path("star", type: "dir")

    """
    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${star_index} \
    --readFilesIn ${fastp_fastq} \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix star/${id}. \
    --twopassMode Basic \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1
    """
}

process align_paired_reads_star {
    label "lrg"
    publishDir "${params.output_directory}/${id}", \
        mode: "copy", overwrite: true

    input:
    tuple val(id), path(fastp_fastqs)
    path(star_index)

    output:
    tuple val(id), path("star", type: "dir")

    """
    STAR \
    --runThreadN ${task.cpus} \
    --genomeDir ${star_index} \
    --readFilesIn ${fastp_fastqs} \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix star/${id}. \
    --twopassMode Basic \
    --quantMode TranscriptomeSAM GeneCounts \
    --outSAMunmapped Within \
    --outFilterMultimapNmax 1
    """
}

process quantify_single_reads_rsem {
    label "lrg"
    publishDir "${params.output_directory}/${id}", \
        mode: "copy", overwrite: true

    input:
    tuple val(id), val(strand), path(aligned_bam)
    path(rsem_reference)

    output:
    tuple val(id), path("rsem", type: "dir")

    """
    case ${strand} in
        *F ) strand=forward;;
        *R ) strand=reverse;;
        *  ) strand=none;;
    esac

    mkdir rsem

    rsem-calculate-expression \
    --estimate-rspd \
    --strandedness \$strand \
    --append-names \
    --num-threads ${task.cpus} \
    --alignments \
    ${aligned_bam} ${rsem_reference}/reference rsem/${id}
    
    awk -v OFS=\$'\\t' '
        NR==1 {
            print "transcript_id",
            "transcript_name",
            "gene_id",
            "gene_name",
            "length",
            "effective_length",
            "expected_count",
            "TPM",
            "FPKM",
            "IsoPct"
        }
        NR!=1 {
            split(\$1,a,"_");
            split(\$2,b,"_");
            print a[1],a[2],b[1],b[2],\$3,\$4,\$5,\$6,\$7,\$8
        }
    ' rsem/${id}.isoforms.results > rsem/${id}.transcripts.results
    
    rm rsem/${id}.isoforms.results
    """
}

process quantify_paired_reads_rsem {
    label "lrg"
    publishDir "${params.output_directory}/${id}", \
        mode: "copy", overwrite: true

    input:
    tuple val(id), val(strand), path(aligned_bam)
    path(rsem_reference)

    output:
    tuple val(id), path("rsem", type: "dir")

    """
    case ${strand} in
        *F ) strand=forward;;
        *R ) strand=reverse;;
        *  ) strand=none;;
    esac

    mkdir rsem

    rsem-calculate-expression \
    --estimate-rspd \
    --strandedness \$strand \
    --append-names \
    --num-threads ${task.cpus} \
    --alignments \
    --paired-end \
    ${aligned_bam} ${rsem_reference}/reference rsem/${id}
    
    awk -v OFS=\$'\\t' '
        NR==1 {
            print "transcript_id",
            "transcript_name",
            "gene_id",
            "gene_name",
            "length",
            "effective_length",
            "expected_count",
            "TPM",
            "FPKM",
            "IsoPct"
        }
        NR!=1 {
            split(\$1,a,"_");
            split(\$2,b,"_");
            print a[1],a[2],b[1],b[2],\$3,\$4,\$5,\$6,\$7,\$8
        }
    ' rsem/${id}.isoforms.results > rsem/${id}.transcripts.results
    
    rm rsem/${id}.isoforms.results
    """
}

process compute_qc_metrics_picard {
    label "lrg"
    publishDir "${params.output_directory}/${id}/qc", \
        mode: "copy", overwrite: true

    input:
    tuple val(id), val(strand), path(star_genome_bam)
    path(ref_flat)
    path(rRNA_interval_list)

    output:
    tuple val(id), path("${id}.rna_metrics")

    """
    case ${strand} in
        *F ) strand=FIRST_READ_TRANSCRIPTION_STRAND;;
        *R ) strand=SECOND_READ_TRANSCRIPTION_STRAND;;
        *  ) strand=NONE;;
    esac

    picard -Xmx16g CollectRnaSeqMetrics \
    --INPUT ${star_genome_bam} \
    --OUTPUT ${id}.rna_metrics \
    --REF_FLAT ${ref_flat} \
    --RIBOSOMAL_INTERVALS ${rRNA_interval_list} \
    --STRAND_SPECIFICITY \$strand \
    --USE_JDK_DEFLATER true \
    --USE_JDK_INFLATER true
    """
}

process count_genes_with_reads {
    label "lrg"

    input:
    tuple val(id), val(strand), path(star_reads_per_gene)

    output:
    tuple val(id), env(n_genes)

    """
    case ${strand} in
        *F ) s=f;;
        *R ) s=r;;
        *  ) s=u;;
    esac

    if [[ \$s == "u" ]]; then
        cut -f1,2 ${star_reads_per_gene} > counts.tsv;
    elif [[ \$s == "f" ]]; then
        cut -f1,3 ${star_reads_per_gene} > counts.tsv;
    elif [[ \$s == "r" ]]; then
        cut -f1,4 ${star_reads_per_gene} > counts.tsv;
    fi

    n_genes=\$(
        awk '
        BEGIN { n=0; }
        NR>=5 { if (\$2 > 0) n+=1; }
        END { print n; }
        ' counts.tsv
    )
    """
}

process merge_qc_metrics {
    label "lrg"
    publishDir "${params.output_directory}/${id}/qc", \
        mode: "copy", overwrite: true

    input:
    tuple val(id), val(n_genes), path(fastp_json), path(star_genome_bam), path(picard_metrics)

    output:
    tuple val(id), path("${id}.qc_metrics.tsv")

    """
    declare -A metrics
    metrics[MAPPED_READS]=\$(samtools view -c -F 4 ${star_genome_bam})
    metrics[N_MAPPED_GENES]=${n_genes}
    metrics[RAW_READS]=\$(cat ${fastp_json} | jq '.summary.before_filtering.total_reads')
    metrics[RAW_BASES]=\$(cat ${fastp_json} | jq '.summary.before_filtering.total_bases')
    metrics[RAW_Q30_RATE]=\$(cat ${fastp_json} | jq '.summary.before_filtering.q30_rate')
    metrics[RAW_READ1_MEAN_LENGTH]=\$(cat ${fastp_json} | jq '.summary.before_filtering.read1_mean_length')
    metrics[RAW_READ2_MEAN_LENGTH]=\$(cat ${fastp_json} | jq '.summary.before_filtering.read2_mean_length')
    metrics[RAW_GC_CONTENT]=\$(cat ${fastp_json} | jq '.summary.before_filtering.gc_content')
    metrics[CLEANED_READS]=\$(cat ${fastp_json} | jq '.summary.after_filtering.total_reads')
    metrics[CLEANED_BASES]=\$(cat ${fastp_json} | jq '.summary.after_filtering.total_bases')
    metrics[CLEANED_Q30_RATE]=\$(cat ${fastp_json} | jq '.summary.after_filtering.q30_rate')
    metrics[CLEANED_READ1_MEAN_LENGTH]=\$(cat ${fastp_json} | jq '.summary.after_filtering.read1_mean_length')
    metrics[CLEANED_READ2_MEAN_LENGTH]=\$(cat ${fastp_json} | jq '.summary.after_filtering.read2_mean_length')
    metrics[CLEANED_GC_CONTENT]=\$(cat ${fastp_json} | jq '.summary.after_filtering.gc_content')

    printf '%s\\t%s\\n' "SAMPLE" ${id} > metrics.tsv
    for key in "\${!metrics[@]}"; do
        printf '%s\\t%s\\n' "\$key" "\${metrics[\$key]}" >> metrics.tsv;
    done

    awk -v FS=\$'\\t' 'NF>3' ${picard_metrics} | \
    cut -f1-27 | \
    awk -v FS=\$'\\t' '
        NR==1 { for (i=1;i<=NF;i++) a[i]=\$i; }
        NR==2 { for (i=1;i<=NF;i++) b[i]=\$i; }
        END { for (i in a) print a[i]"\\t"b[i] }
    ' > picard.tsv

    cat metrics.tsv picard.tsv > ${id}.qc_metrics.tsv
    """
}

process collect_salmon_quants {
    label "lrg"
    publishDir "${params.output_directory}/summary", \
        mode: "copy", overwrite: true

    input:
    path(gtf)
    path(salmon_quants)

    output:
    path("${params.dataset_name}.salmon.*.*.tximport.tsv")

    """
    collect_quants.R \
    ${params.dataset_name} \
    "salmon" \
    ${gtf} \
    ${salmon_quants} 
    """
}

process collect_rsem_quants {
    label "lrg"
    publishDir "${params.output_directory}/summary", \
        mode: "copy", overwrite: true

    input:
    path(gtf)
    path(rsem_quants)

    output:
    path("${params.dataset_name}.rsem.*.*.tximport.tsv")

    """
    collect_quants.R \
    ${params.dataset_name} \
    "rsem" \
    ${gtf} \
    ${rsem_quants}
    """
}

process collect_qc_metrics {
    label "lrg"
    publishDir "${params.output_directory}/summary", \
        mode: "copy", overwrite: true

    input:
    path(qc_metrics)

    output:
    path("${params.dataset_name}.qc_metrics.tsv")

    """
    awk -v FS=\$'\\t' '
            NR==FNR {
                header=header "\\t" \$1;
            }
            FNR==1 {
                gsub(".qc_metrics.tsv", "", FILENAME);
                a[FILENAME]=\$2;
            }
            FNR>1 { 
                a[FILENAME]=a[FILENAME] "\\t" \$2;
            }
            END {
                print header;
                for (sample in a) {
                    print a[sample];
                }
            }
    ' *.qc_metrics.tsv > ${params.dataset_name}.qc_metrics.tsv
    """
}

process generate_multiqc_report {
    label "lrg"
    publishDir "${params.output_directory}/summary", \
        mode: "copy", overwrite: true

    input:
    path(multiqc_config)
    path(pipeline_logs)

    output:
    path("${params.dataset_name}.multiqc_report.html")

    """
    multiqc \
    --config ${multiqc_config} \
    --filename ${params.dataset_name}.multiqc_report.html \
    ${pipeline_logs}
    """
}
