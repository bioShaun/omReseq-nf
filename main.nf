#!/usr/bin/env nextflow

def helpMessage() {
    log.info """

    Usage:

    nextflow run omReseq-nf/ --reads '*_R{1,2}.fastq.gz' --genome GRCh37 -profile uppmax

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Other options:
      --outdir                      The output directory where the results will be saved

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

 