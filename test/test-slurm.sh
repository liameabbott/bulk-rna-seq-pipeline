#!/bin/bash
export AWS_PROFILE=liam
nextflow run $(pwd)/../main.nf -profile slurm,conda,test
