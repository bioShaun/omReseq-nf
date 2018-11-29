#!/usr/bin/env nextflow

def helpMessage() {
    log.info """

    Usage:

    =======================================================
    nextflow run omReseq-nf --reads '*_R{1,2}.fastq.gz'
    =======================================================

    Test arguments:
      --test                        The project is for test, default is false
      --test_data_size              Test data size, default is 100000
      --test_data_dir               Test data directory

    References If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Mandatory arguments:
      --reads_dir                   Path to input data (must be surrounded with quotes)

    Other options:
      --outdir                      The output directory where the results will be saved

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

 // Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// default parameters
resds_pattern = '*.R{1,2}.clean.fastq.gz'
// resds_pattern = '*_{1,2}.fq'
params.skip_qc = false
params.skip_fastqc = false
params.test = true
params.test_data_size = 1000 * 1000
params.test_data_dir = false
params.fasta = false

// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads_dir}/${resds_pattern}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${resds_pattern} in ${params.reads_dir}!" }
    .into { raw_reads }


process fetch_fastq {
    tag "Fetch fq on $name"

    publishDir "${params.outdir}/fq_dir"

    input:
    set name, file(reads) from raw_reads

    output:
    file "*fq.gz" into fastqc_fq_files, bwa_fq_files

    script:
    if (params.test) {
        """
        seqtk sample -s100 ${reads[0]} ${params.test_data_size} | gzip > ${name}.R1.fq.gz
        seqtk sample -s100 ${reads[1]} ${params.test_data_size} | gzip > ${name}.R2.fq.gz
        """
    } else {
        """
        ln -s ${reads[0]} ${name}.R1.fq.gz
        ln -s ${reads[1]} ${name}.R2.fq.gz
        """
    }

}



/*
 * STEP 0 - FastQC
 */
process fastqc {
    tag "FASTQC on $name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    !params.skip_qc && !params.skip_fastqc

    input:
    file reads from fastqc_fq_files

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    cpus = 8

    script:
    name = reads[0].toString() - ~/.R1.fq.gz$/
    """
    fastqc -q ${reads}
    """
}


/*
* BWA Mapping
*/

process bwa_mapping {
    tag "BWA Mapping on ${sample_name}"
    publishDir "${params.outdir}/mapping/${sample_name}"

    input:
    file reads from bwa_fq_files

    output:
    file "${sample_name}.sort.bam" into samtools_stats_bam, picard_bam

    cpus = 16

    script:    
    sample_name = reads[0].toString() - ~/.R1.fq.gz$/
    """
    bwa mem -M -a \
	    -R \"@RG\tID:${sample_name}\tSM:${sample_name}\tLB:${sample_name}\tPI:350\tPL:Illumina\tCN:TCuni\" \
	    -t ${task.cpus} \
	    -K 10000000 \
	    ${params.fasta} \
	    ${reads[0]} \
	    ${reads[0]} \
	    | \
	samtools view -O bam \
	    --threads ${task.cpus} \
	    -o ${sample_name}.bam
	
    samtools sort -m 2400M --threads ${task.cpus} \
	    -o ${sample_name}.sort.bam \
	    ${sample_name}.bam
    """
}
