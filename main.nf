#!/usr/bin/env nextflow

// TODO
// 3. test fq generate process
// 4. exon/cds bed generate

def helpMessage() {
    log.info """

    Usage:

    =======================================================
    nextflow run omReseq-nf --reads '*_R{1,2}.fastq.gz'
    =======================================================

    References If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --bwa_index                   Path to reference bwa index
      --exon_bed                    Path to reference exon bed file
      --cds_bed                     Path to reference cds bed file

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

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
params.skip_qc = true
params.skip_fastqc = false
params.fasta = false
params.bwa_index = false
params.reads = false
params.cds_bed = false
params.exon_bed = false

// reference files
cds_bed_file = file(params.cds_bed)
exon_bed_file = file(params.exon_bed)
bwa_index_file = file(params.bwa_index)
genome_fa = file(params.fasta)


// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${resds_pattern} in ${params.reads_dir}!" }
    .into { fastqc_fq_files;  bwa_fq_files }

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
    set name, file(reads) from fastqc_fq_files

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    cpus = 2

    script:
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
    set sample_name, file(reads) from bwa_fq_files
    file bwa_index_file from bwa_index_file
    file genome_fa from genome_fa

    output: 
    file "${sample_name}.sort.bam" into samtools_stats_bam, picard_bam

    cpus = 16

    script:    
    """
    bwa mem -M -a \
	    -R \"@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPI:350\\tPL:Illumina\\tCN:TCuni\" \
	    -t ${task.cpus} \
	    -K 10000000 \
	    ${bwa_index_file}/${genome_fa.getName()} \
	    ${reads[0]} \
	    ${reads[1]} \
	    | \
	samtools view -O bam \
	    --threads ${task.cpus} \
	    -o ${sample_name}.bam
	
    samtools sort -m 2400M --threads ${task.cpus} \
	    -o ${sample_name}.sort.bam \
	    ${sample_name}.bam
    """
}


process reads_cov_stats {

    tag "SAMTOOLS stats on ${sample_name}"

    publishDir "${params.outdir}/mapping/${sample_name}"

    input:
    file samtools_stats_bam from samtools_stats_bam
    file cds_bed_file from cds_bed_file
    file exon_bed_file from exon_bed_file

    output:
    file "${sample_name}.*.stat" into reads_cov_results

    cpus = 8

    script:
    sample_name = samtools_stats_bam.baseName - '.sort'
    """
    samtools stats \
	    --threads ${task.cpus} \
	    --target-regions  ${cds_bed_file} \
	    ${samtools_stats_bam} \
	    > ${sample_name}.cds.stat

    samtools stats \
	    --threads ${task.cpus} \
	    --target-regions  ${exon_bed_file} \
	    ${samtools_stats_bam} \
	    > ${sample_name}.exon.stat

    samtools stats \
	    --threads 20 \
	    ${samtools_stats_bam} \
	    > ${sample_name}.genome.stat 
    """
}