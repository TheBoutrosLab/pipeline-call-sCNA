#!/usr/bin/env python3
"""
- This script is used to run SigProfilerExtractor to extract mutational signatures from a mutational
catalog matrix.
- The input is a tab delimited file containing the mutational catalog matrix
- The output is a folder containing the extracted mutational signatures
"""

import argparse
from SigProfilerExtractor import sigpro as sig

parser = argparse.ArgumentParser(
    description='Run SigProfilerExtractor to extract mutational signatures \
        from a mutational catalog matrix'
    )
parser.add_argument(
    '--cna-caller',
    type=str,
    help='CNA caller used to generate input'
    )
parser.add_argument(
    '--sample',
    type=str,
    help='Sample name'
    )
parser.add_argument(
    '--input',
    type=str,
    help='Path to the input file containing the mutational catalog matrix'
    )
parser.add_argument(
    '--output-dir',
    type=str,
    default='.',
    help='Base output directory (default: current directory)'
    )
parser.add_argument(
    '--reference-genome',
    type=str,
    default='GRCh38',
    help='Reference genome. This parameter is applicable only if the input_type is "vcf"'
    )
parser.add_argument(
    '--opportunity-genome',
    type=str,
    help='The build or version of the reference genome for the reference signatures'
    )
parser.add_argument(
    '--context-type',
    type=str,
    default='default',
    help='Mutation Context type for signature extraction'
    )
parser.add_argument(
    '--exome',
    action='store_true',
    help='Enable exome mode'
    )
parser.add_argument(
    '--min-signatures',
    type=int,
    default=1,
    help='Minimum number of signatures to be extracted'
    )
parser.add_argument(
    '--max-signatures',
    type=int,
    default=25,
    help='Maximum number of signatures to be extracted'
    )
parser.add_argument(
    '--nmf-replicates',
    type=int,
    default=100,
    help='Number iteration to be performed to extract each number signature'
    )
parser.add_argument(
    '--resample',
    action='store_true',
    help='Add poisson noise to samples by resampling'
    )
parser.add_argument(
    '--seeds',
    type=str,
    default='random',
    help='Ensures reproducible NMF replicate resamples'
    )
parser.add_argument(
    '--matrix-normalization',
    type=str,
    default='gmm',
    choices=['gmm','log2', 'custom', 'none'],
    help='Method of normalizing the genome matrix before it is analyzed by NMF'
    )
parser.add_argument(
    '--nmf-init',
    type=str,
    default='random',
    choices=['random', 'nndsvd', 'nndsvda', 'nndsvdar', 'nndsvd_min'],
    help='The initialization algorithm for W and H matrix of NMF'
    )
parser.add_argument(
    '--precision',
    type=str,
    default='single',
    choices=['single','double'],
    help='Values should be single or double'
    )
parser.add_argument(
    '--min-nmf-iterations',
    type=int,
    default=10000,
    help='Value defines the minimum number of iterations to be completed before NMF converges'
    )
parser.add_argument(
    '--max-nmf-iterations',
    type=int,
    default=1000000,
    help='Value defines the maximum number of iterations to be completed before NMF converges'
    )
parser.add_argument(
    '--nmf-test-conv',
    type=int,
    default=10000,
    help='Value defines the number number of iterations to done between checking next convergence'
    )
parser.add_argument(
    '--nmf-tolerance',
    type=float,
    default=1e-15,
    help='Value defines the tolerance to achieve to converge'
    )
parser.add_argument(
    '--cpu',
    type=int,
    default=-1,
    help='The number of processors to be used to extract the signatures (default: -1 for all)'
    )
parser.add_argument(
    '--gpu',
    action='store_true',
    help='Defines if the GPU resource will used if available'
    )
parser.add_argument(
    '--batch-size',
    type=int,
    default=1,
    help='Defines the number of NMF replicates to be performed \
        by each CPU during the parallel processing'
    )
parser.add_argument(
    '--stability',
    type=float,
    default=0.8,
    help='The cutoff thresh-hold of the average stability. \
        Solutions with average stabilities below this thresh-hold will not be considered.'
    )
parser.add_argument(
    '--min-stability',
    type=float,
    default=0.2,
    help='The cutoff thresh-hold of the minimum stability. \
        Solutions with minimum stabilities below this thresh-hold will not be considered.'
    )
parser.add_argument(
    '--combined-stability',
    type=float,
    default=1.0,
    help='The cutoff thresh-hold of the combined stability (sum of average and minimum stability).\
        Solutions with combined stabilities below this thresh-hold will not be considered.'
    )
parser.add_argument(
    '--allow-stability-drop',
    action='store_true',
    help='Defines if solutions with a drop in stability with respect to \
        the highest stable number of signatures will be considered'
    )
parser.add_argument(
    '--cosmic-version',
    type=float,
    default='3.4',
    help='Defines the version of the COSMIC reference signatures'
    )
parser.add_argument(
    '--make-decomposition-plots',
    action='store_true',
    help='Generate de novo to COSMIC signature decomposition plots as part of the results'
    )
parser.add_argument(
    '--collapse-to-SBS96',
    action='store_true',
    help='If True, SBS288 and SBS1536 de novo signatures will be mapped \
        to SBS96 reference signatures. \
        If False, those will be mapped to reference signatures of the same context.'
    )
parser.add_argument(
    '--get-all-signature-matrices',
    action='store_true',
    help='Write to output Ws and Hs from all the NMF iterations'
    )
parser.add_argument(
    '--export-probabilities',
    action='store_true',
    help='Create the probability matrix'
    )
parser.add_argument(
    '--volume',
    type=str,
    default='.',
    help='Path to the volume for writing and loading reference genomes,\
        plotting templates, and COSMIC signature plots'
    )

args = parser.parse_args()

def main_function():
    """
    Function to run SigProfilerExtractor
    """

    sig.sigProfilerExtractor(
    f'seg:{args.cna_caller}',
    args.output_dir,
    args.input,
    reference_genome=args.reference_genome,
    opportunity_genome=args.opportunity_genome,
    context_type=args.context_type,
    exome=args.exome,
    minimum_signatures=args.min_signatures,
    maximum_signatures=args.max_signatures,
    nmf_replicates=args.nmf_replicates,
    resample=args.resample,
    batch_size=args.batch_size,
    cpu=args.cpu,
    gpu=args.gpu,
    nmf_init=args.nmf_init,
    precision=args.precision,
    matrix_normalization=args.matrix_normalization,
    seeds=args.seeds,
    min_nmf_iterations=args.min_nmf_iterations,
    max_nmf_iterations=args.max_nmf_iterations,
    nmf_test_conv=args.nmf_test_conv,
    nmf_tolerance=args.nmf_tolerance,
    stability=args.stability,
    min_stability=args.min_stability,
    combined_stability=args.combined_stability,
    allow_stability_drop=args.allow_stability_drop,
    cosmic_version=args.cosmic_version,
    make_decomposition_plots=args.make_decomposition_plots,
    collapse_to_SBS96=args.collapse_to_SBS96,
    get_all_signature_matrices=args.get_all_signature_matrices,
    export_probabilities=args.export_probabilities,
    volume=args.output_dir
    )

if __name__=="__main__":
    main_function()
