/* MIT License

 * Copyright (c) 2021, Liam Abbott.

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// Enable DSL2 syntax.
nextflow.enable.dsl = 2

// List of required parameters.
def required_params = [
    "multiqc_config",
    "dataset_name",
    "fastq_directory",
    "output_directory",
    "reference_fasta",
    "reference_gtf",
    "salmon_index",
    "star_index",
    "rsem_reference",
    "reference_refFlat",
    "rRNA_interval_list"
]

// Check that all required params are present.
for (param in required_params) {
    if ( !params[param] ) {
      error "Missing ${param} parameter."
    }
}

// Import pipeline modules.
include {
    clean_single_reads_fastp;
    clean_paired_reads_fastp;
    quantify_single_reads_salmon;
    quantify_paired_reads_salmon;
    align_single_reads_star;
    align_paired_reads_star;
    quantify_single_reads_rsem;
    quantify_paired_reads_rsem;
    compute_qc_metrics_picard;
    count_genes_with_reads;
    merge_qc_metrics;
    collect_salmon_quants;
    collect_rsem_quants;
    collect_qc_metrics;
    generate_multiqc_report
} from './modules.nf'

// Main pipeline logic.
workflow {

  // Load single- and paired-end reads from directory.
  reads = Channel
    .fromFilePairs("${params.fastq_directory}/*.fastq.gz", flat: true, size: -1)
    .map{ it -> [it[0].split("_")[0], it[1]] }
    .groupTuple()
    .branch {
      single_end: it[1].size() == 1
      paired_end: it[1].size() >= 2
    }

  // Logic to split fastp output
  def fastp_criteria = multiMapCriteria { it ->
    fastqs: [it[0], it[1].findAll{ f -> f =~ /.fastp.fastq.gz/ }]
    json: [it[0], it[1].find{ f -> f =~ /.fastp.json/ }]
  }

  // Clean raw reads.
  fastp_single = clean_single_reads_fastp(reads.single_end)
    .multiMap(fastp_criteria) 
  fastp_paired = clean_paired_reads_fastp(reads.paired_end)
    .multiMap(fastp_criteria)
  fastp_jsons = fastp_single.json.mix(fastp_paired.json)

  // Load salmon index.
  salmon_index = Channel
    .fromPath(params.salmon_index, type: "dir", checkIfExists: true)
    .collect()

  // Load reference GTF annotation file.
  reference_gtf = Channel
    .fromPath(params.reference_gtf, checkIfExists: true)
    .collect()

  // Logic to split salmon output.
  def salmon_criteria = multiMapCriteria { it ->
    strand: [it[0], it[1]]
    transcript_quant: [it[0], it[2] + "/${it[0]}.quant.sf"]
    gene_quant: [it[0], it[2] + "/${it[0]}.quant.genes.sf"]
    directory: [it[0], it[2]]
  }
  
  // Quantification using salmon.
  salmon_single = quantify_single_reads_salmon(fastp_single.fastqs, salmon_index, reference_gtf)
    .multiMap(salmon_criteria)
  salmon_paired = quantify_paired_reads_salmon(fastp_paired.fastqs, salmon_index, reference_gtf)
    .multiMap(salmon_criteria)

  // Merge single- and paired-end salmon output.
  salmon_strands = salmon_single.strand.mix(salmon_paired.strand)
  salmon_transcript_quants = salmon_single.transcript_quant.mix(salmon_paired.transcript_quant)
  salmon_gene_quants = salmon_single.gene_quant.mix(salmon_paired.gene_quant)
  salmon_directories = salmon_single.directory.mix(salmon_paired.directory)

  // Load STAR index.
  star_index = Channel
    .fromPath(params.star_index, type: "dir", checkIfExists: true)
    .collect()

  // Logic to split STAR output.
  def star_criteria = multiMapCriteria { it ->
    genome_bam: [it[0], it[1] + "/${it[0]}.Aligned.sortedByCoord.out.bam"]
    transcriptome_bam: [it[0], it[1] + "/${it[0]}.Aligned.toTranscriptome.out.bam"]
    reads_per_gene: [it[0], it[1] + "/${it[0]}.ReadsPerGene.out.tab"]
    log: [it[0], it[1] + "/${it[0]}.Log.final.out"]
  }

  // Align reads.
  star_single = align_single_reads_star(fastp_single.fastqs, star_index)
    .multiMap(star_criteria)
  star_paired = align_paired_reads_star(fastp_paired.fastqs, star_index)
    .multiMap(star_criteria)
  
  // Merge single- and paired-end STAR output.
  star_genome_bams = star_single.genome_bam.mix(star_paired.genome_bam)
  star_reads_per_gene = star_single.reads_per_gene.mix(star_paired.reads_per_gene)  
  star_logs = star_single.log.mix(star_paired.log)

  // Load RSEM reference index.
  rsem_reference = Channel
    .fromPath(params.rsem_reference, type: "dir", checkIfExists: true)
    .collect()

  // Logic to split RSEM output.
  def rsem_criteria = multiMapCriteria { it -> 
    transcript_quant: [it[0], it[1] + "/${it[0]}.transcripts.results"]    
    gene_quant: [it[0], it[1] + "/${it[0]}.genes.results"]
    log: [it[0], it[1] + "/${it[0]}.stat/${it[0]}.cnt"]
  }

  // Quantify aligned reads.
  rsem_single = quantify_single_reads_rsem(
    salmon_single.strand.join(star_single.transcriptome_bam), rsem_reference).multiMap(rsem_criteria)
  rsem_paired = quantify_paired_reads_rsem(
    salmon_paired.strand.join(star_paired.transcriptome_bam), rsem_reference).multiMap(rsem_criteria)

  // Merge single- and paired-end RSEM output.
  rsem_transcript_quants = rsem_single.transcript_quant.mix(rsem_paired.transcript_quant)
  rsem_gene_quants = rsem_single.gene_quant.mix(rsem_paired.gene_quant)
  rsem_logs = rsem_single.log.mix(rsem_paired.log) 

  // Load refFlat annotation file.
  ref_flat = Channel
    .fromPath(params.reference_refFlat, checkIfExists: true)
    .collect()

  // Load rRNA interval list.
  rRNA_interval_list = Channel
    .fromPath(params.rRNA_interval_list, checkIfExists: true)
    .collect()

  // Compute QC metrics.
  picard_metrics = compute_qc_metrics_picard(
    salmon_strands.join(star_genome_bams), ref_flat, rRNA_interval_list)
  gene_counts = count_genes_with_reads(
    salmon_strands.join(star_reads_per_gene))
  qc_metrics = merge_qc_metrics(
    gene_counts.join(fastp_jsons).join(star_genome_bams).join(picard_metrics))

  // Collect quants and QC metrics for all samples.
  collect_salmon_quants(
    reference_gtf, salmon_transcript_quants.map{ it -> it[1] }.collect())
  collect_rsem_quants(
    reference_gtf, rsem_transcript_quants.map { it -> it[1] }.collect())  
  collect_qc_metrics(
    qc_metrics.map{ it -> it[1] }.collect())
  
  // Load multiqc configuration file.
  multiqc_config = Channel.fromPath("${params.multiqc_config}")

  // Generate multiqc report.
  generate_multiqc_report(
    multiqc_config,    
    fastp_jsons.map{ it -> it[1] } 
      .mix(salmon_directories.map{ it -> it[1] })
      .mix(star_logs.map{ it -> it[1] })
      .mix(star_reads_per_gene.map{ it -> it[1] })
      .mix(rsem_logs.map{ it -> it[1] })
      .mix(picard_metrics.map{ it -> it[1] }).collect())

}

