#!/usr/bin/env nextflow

// TODO
// No cds bed file?


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
      --split_bed                   Path to split bed directory
      --known_vcf                   Path to known vcf file

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
params.known_vcf = false
params.split_bed = false

// reference files
cds_bed_file = file(params.cds_bed)
exon_bed_file = file(params.exon_bed)
bwa_index_file = file(params.bwa_index)
genome_fa = file(params.fasta)
genome_path = genome_fa.getParent()
genome_fai = file("${params.fasta}.fai")
genome_dict = file("${genome_path}/${genome_fa.baseName}.dict")
known_vcf = file(params.known_vcf)
known_vcf_index = file("${params.known_vcf}.tbi")
split_bed_files = Channel.fromPath("${params.split_bed}/*bed")

// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${resds_pattern} in ${params.reads_dir}!" }
    .into { fastqc_fq_files;  bwa_fq_files }

/*
 * FastQC
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
    file "${sample_name}.bam" into unsort_bam

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
    """
}


/*
* Sort bam
*/
process sort_bam {
    tag "Sort Bam on ${sample_name}"

    publishDir "${params.outdir}/mapping/${sample_name}"

    input:
    file bam from unsort_bam

    output: 
    file "${sample_name}.sort.bam" into samtools_stats_bam, to_rmdup_bam

    cpus = 8

    script:
    sample_name = bam.baseName
    """	
    samtools fixmate -m ${bam} ${sample_name}.fixmate.bam

    samtools sort -m 2400M --threads ${task.cpus} \
	    -o ${sample_name}.sort.bam \
	    ${sample_name}.fixmate.bam
    """
}

/*
* Reads coverage stats
*/
process reads_cov_stats {

    tag "SAMTOOLS Stats on ${sample_name}"

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
	    --threads ${task.cpus} \
	    ${samtools_stats_bam} \
	    > ${sample_name}.genome.stat 
    """
}


/*
* Reads mark duplication and recalibration
*/

process bam_remove_duplicate {
    tag "REMOVE DUP on ${sample_name}"

    publishDir "${params.outdir}/mapping/${sample_name}"

    input:
    file bam from to_rmdup_bam
  
    output:
    file "${sample_name}.rmdup.bam" into br_rmdup_bam
  
    cpus = 8

    script:
    sample_name = bam.baseName - '.sort'
    """
    samtools markdup -r -s \
	    --threads ${task.cpus} \
	    ${bam} \
	    ${sample_name}.rmdup.bam
    """
}


/*
* bam BaseRecalibrator
*/
process bam_BaseRecalibrator {
    tag "BaseRecalibrator on ${sample_name}"

    publishDir "${params.outdir}/mapping/${sample_name}"

    input:
    file bam from br_rmdup_bam
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict
    file known_vcf from known_vcf
    file known_vcf_index from known_vcf_index
  
    output:
    file "${sample_name}.recal.table" into recal_table
    file bam into bqsr_rmdup_bam
  
    cpus = 8

    script:
    sample_name = bam.baseName - '.rmdup'
    """
    gatk BaseRecalibrator \
        --reference ${refer} \
        --input ${sample_name}.rmdup.bam \
        --output ${sample_name}.recal.table \
        --known-sites ${known_vcf}
    """
}

/*
* bam ApplyBQSR
*/
process bam_ApplyBQSR {
    tag "ApplyBQSR on ${sample_name}"

    publishDir "${params.outdir}/mapping/${sample_name}"

    input:
    file bam from bqsr_rmdup_bam
    file recal_table_file from recal_table
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict
  
    output:
    file "${sample_name}.bqsr.bam" into bqsr_bam
  
    cpus = 8

    script:
    sample_name = bam.baseName - '.rmdup'
    """       
    gatk ApplyBQSR \
        --bqsr-recal-file ${recal_table_file} \
        --input ${bam} \
        --output ${sample_name}.bqsr.bam \
        --read-filter AmbiguousBaseReadFilter \
        --read-filter MappingQualityReadFilter \
        --read-filter NonZeroReferenceLengthAlignmentReadFilter \
        --read-filter ProperlyPairedReadFilter \
        --minimum-mapping-quality 30 
    """
}

/*
*  GATK HaplotypeCaller
*/

process gatk_HaplotypeCaller {
    tag "GATK HaplotypeCaller on ${sample_name} - ${chr_name}"

    publishDir "${params.outdir}/gvcf/${sample_name}"

    input:
    file bam from bqsr_bam
    each file(bed) from split_bed_files
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    

    output:
    file "${sample_name}.${chr_name}.hc.g.vcf.gz" into sample_gvcf
    
    cpus = 8

    script:
    sample_name = bam.baseName - '.bqsr'
    chr_name = bed.baseName
    """
    gatk HaplotypeCaller  \
        --input ${bam} \
        --output ${sample_name}.${chr_name}.hc.g.vcf.gz \\
        --reference ${refer} \\
        --intervals ${bed} \\
        --emit-ref-confidence GVCF
    """
}