#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(tximport)
library(GenomicFeatures)

args = commandArgs(trailingOnly=TRUE)
dataset_name = args[1]
quant_type = args[2]
gtf = args[3]
quants = args[4:length(args)]

txdb = makeTxDbFromGFF(gtf)
k = keys(txdb,  keytype='TXNAME')
tx2gene = AnnotationDbi::select(
  txdb, k, 'GENEID', 'TXNAME')

if (quant_type == 'salmon') {
    filename = '(.*)\\.quant\\.sf'
    gene_id_col = NULL
    tx_id_col = 'Name'
    abundance_col = 'TPM'
    counts_col = 'NumReads' 
    length_col = 'EffectiveLength'
} else if (quant_type == 'rsem') {
    filename = '(.*)\\.transcripts\\.results'
    gene_id_col = 'gene_id'
    tx_id_col = 'transcript_id'
    abundance_col = 'TPM'
    counts_col = 'expected_count'
    length_col = 'effective_length'
}

importer = function(x) readr::read_tsv(x, progress=FALSE)
names(quants) = stringr::str_match(quants, filename)[,2]

txi.tx = tximport(
    quants, type='none', txIn=TRUE, txOut=TRUE, importer=importer,
    geneIdCol=gene_id_col, txIdCol=tx_id_col, abundanceCol=abundance_col,
    countsCol=counts_col, lengthCol=length_col)

txi.gene = summarizeToGene(
  txi.tx, tx2gene, ignoreTxVersion=TRUE)

extract_field_df = function(field, name, feature_type) {
  df = rownames_to_column(data.frame(field))
  id_col = paste(
    substr(feature_type, 1, nchar(feature_type)-1),
    'id', sep='_')
  n_samples = dim(df)[2]
  colnames(df) = c(id_col, colnames(df)[2:n_samples])
  df = df[,c(1, order(colnames(df)[2:n_samples])+1)]
  write.table(
    df, row.names=FALSE, sep='\t', quote=FALSE,
    file=paste0(
      dataset_name, '.', quant_type, '.',
      feature_type, '.', name, '.tximport.tsv'))
}

extract_field_df(txi.tx$length, 'effective_length', 'transcripts')
extract_field_df(txi.tx$abundance, 'TPM', 'transcripts')
extract_field_df(txi.tx$counts, 'read_count', 'transcripts')
extract_field_df(txi.gene$length, 'effective_length', 'genes')
extract_field_df(txi.gene$abundance, 'TPM', 'genes')
extract_field_df(txi.gene$counts, 'read_count', 'genes')

