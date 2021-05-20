#!/usr/bin/env python

import os
import sys
import argparse
from subprocess import Popen, PIPE, STDOUT

def main(args_dict):
    
    # Run pipeline with defined parameters.
    nextflow_pipeline = Popen([
        'nextflow',
        'run',
        '/shared/bulk-rna-seq-processing-pipeline/main.nf',
        '-latest',
        '-profile', args_dict['nxf_profile'],
        '--dataset-name', args_dict['dataset_name'],
        '--reference_database', args_dict['reference_database'],
        '--fastq_directory', args_dict['fastq_directory'],
        '--output_directory', args_dict['output_directory'],
        '--source_database', args_dict['source_database'],
        '--source_database_release', args_dict['source_database_release'],
        '--genome_species', args_dict['genome_species'],
        '--genome_assembly', args_dict['genome_assembly'],
        '--rRNA_interval_list', args_dict['rRNA_interval_list']
    ], stdout=PIPE, stderr=STDOUT)

    # Print pipeline stdout and stderr streams.
    for line in iter(nextflow_pipeline.stdout.readline, b''):
        sys.stdout.write(line.decode(sys.stdout.encoding))

if __name__ == '__main__':

    # Parse command-line arguments.
    parser = argparse.ArgumentParser(
        description='Process raw bulk RNA-seq FASTQs into gene/transcript counts.')
    parser.add_argument(
        '--public-or-private', '-p', choices=['public', 'private'],
        help='Dataset type.')
    parser.add_argument(
        '--project-name', '-n', 
        help='Name of project (e.g. "poi2pod-biomarker").')
    parser.add_argument(
        '--dataset-name', '-t', 
        help='Dataset name (such as SRA bioproject ID, e.g. "PRJNA315611").')
    parser.add_argument(
        '--source-database', '-s', choices=['ensembl', 'gencode', 'refseq'],
        help='Reference database source to use.')
    parser.add_argument(
        '--source-database-release', '-r',
        help='Release/version of reference database source (e.g. "102" for "ensembl").')
    parser.add_argument(
        '--genome-species', '-c',
        help='Species reference genome to use in processing the data.')
    parser.add_argument(
        '--genome-assembly', '-a', 
        help='Assembly version of the reference genome to use.')
    parser.add_argument(
        '--rRNA-interval-list', '-i', default=None,
        help='Path to a specific rRNA interval list to use in analysis.')
    parser.add_argument(
        '--reference-database', '-d',
        default='s3://gates-mri-bioinformatics/reference-database',
        help='Path to the reference database.')
    parser.add_argument(
        '--root-directory', default='s3://gates-mri-bioinformatics/',
        help='Root directory or S3 bucket.')
    parser.add_argument(
        '--nxf-profile', '-x', default='bio', choices=['bio', 'sge'],
        help='Nextflow execution profile.')
    parser.add_argument(
        '--dev-mode', '-v', action='store_true',
        help='Run pipeline with example dev dataset.')
    args = parser.parse_args()

    # Convert arguments into dictionary.
    args_dict = vars(args)

    # Default reference database.
    if not args_dict['reference_database']:
        args_dict['reference_database'] = os.path.join(
            args_dict['root_directory'], 'reference-database')

    # Use specific example dataset if run in dev mode.
    if args_dict['dev_mode']:
        args_dict['dataset_name'] = 'dev'
        args_dict['fastq_directory'] = os.path.join(
            args_dict['root_directory'],
            'PlayGround',
            'liam',
            'bulk-rna-seq-processing-pipeline-dev',
            'sample-reads'
        )
        args_dict['output_directory'] = os.path.join(
            args_dict['root_directory'],
            'PlayGround',
            'liam',
            'bulk-rna-seq-processing-pipeline-dev',
            'processed'
        )
        args_dict['source_database'] = 'ensembl'
        args_dict['source_database_release'] = '102'
        args_dict['genome_species'] = 'homo_sapiens'
        args_dict['genome_assembly'] = 'GRCh38'

    # If production mode.
    else:

        # All provided arguments.
        arg_keys = args_dict.keys()

        # Make sure required arguments are present.
        req_args = [
            'public_or_private',
            'project_name',
            'dataset_name',
            'source_database',
            'source_database_release',
            'genome_species',
            'genome_assembly'
        ]
        for arg in req_args:
            if arg not in arg_keys:
                raise ValueError(f'"{arg}" argument must be provided.')

        # Derived parameters.
        args_dict['fastq_directory'] = os.path.join(
            args_dict['root_directory'],
            f'{args_dict["public_or_private"]}-data',
            args_dict['project_name'],
            args_dict['dataset_name'],
            'fastqs'
        )
        args_dict['output_directory'] = os.path.join(
            args_dict['root_directory'],
            f'{args_dict["public_or_private"]}-data',
            args_dict['dataset_name'],
            args_dict['bioproject_id'],
            'processed'
        )

    # Default rRNA interval list.
    if not args_dict['rRNA_interval_list']:
        args_dict['rRNA_interval_list'] = os.path.join(
            args_dict['reference_database'],
            args_dict['source_database'],
            f'release-{args_dict["source_database_release"]}',
            args_dict['genome_species'],
            args_dict['genome_assembly'],
            'interval-list',
            'reference.genes.rRNA.interval_list'
        )

    # Pass command-line arguments to 
    main(args_dict)
