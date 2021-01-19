#!/usr/bin/env Rscript

library(dplyr)
library(tibble)
library(tximport)
library(GenomicFeatures)
library(readr)

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
  names(quants) = stringr::str_match(
    quants, '(.*)\\.quant\\.sf')[,2]
} else if (quant_type == 'rsem') {
  names(quants) = stringr::str_match(
    quants, '(.*)\\.transcripts.results')[,2]
}

txi.tx = tximport(
  quants, type=quant_type, txIn=TRUE, txOut=TRUE)

txi.gene = summarizeToGene(
  txi.tx, tx2gene, ignoreTxVersion=TRUE)

extract_field_df = function(field, name, feature_type) {
  df = rownames_to_column(data.frame(field))
  colnames(df) = c('transcript_id', colnames(df)[2:dim(df)[2]])
  write.table(
    df, row.names=FALSE, sep='\t', quote=FALSE,
    file=paste0(dataset_name, '.', quant_type, '.', feature_type, '.', name, '.tximport.tsv'))
}

extract_field_df(txi.tx$length, 'effective_length', 'transcripts')
extract_field_df(txi.tx$abundance, 'TPM', 'transcripts')
extract_field_df(txi.tx$counts, 'read_count', 'transcripts')
extract_field_df(txi.gene$length, 'effective_length', 'genes')
extract_field_df(txi.gene$abundance, 'TPM', 'genes')
extract_field_df(txi.gene$counts, 'read_count', 'genes')

