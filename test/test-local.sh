#!/bin/bash
export AWS_PROFILE=liam
nextflow run ../main.nf -profile bio-dev,conda,test
