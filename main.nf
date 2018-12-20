#!/usr/bin/env nextflow

// TODO
// No cds bed file?


def helpMessage() {
    log.info """

    Usage:

    References If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --bwa_index                   Path to reference bwa index
      --exon_bed                    Path to reference exon bed file
      --cds_bed                     Path to reference cds bed file
      --split_bed              Path to split bed directory
      --known_vcf                   Path to known vcf file

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Other options:
      --snpEff_db
      --quality                     
      --depth
      --exome                       Is the project a exom sequencing project
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
// 30 for real data
params.quality = 30
// 5 for read data
params.depth = 5
// snpEff
params.snpEff = '/public/software/snpEff/snpEffv4.3T/'
params.snpEff_db = false
//
params.exome = false

// reference files
cds_bed_file = file(params.cds_bed)
exon_bed_file = file(params.exon_bed)
bwa_index_file = file(params.bwa_index)
genome_fa = file(params.fasta)
genome_path = genome_fa.getParent()
genome_fai = file("${params.fasta}.fai")
genome_dict = file("${genome_path}/${genome_fa.baseName}.dict")
if (params.known_vcf) {
    known_vcf = file(params.known_vcf)
    known_vcf_index = file("${params.known_vcf}.tbi")
} else {
    known_vcf = false
    known_vcf_index = false
}

if (params.exome) {
    split_bed_dir  = "${params.split_bed}/exon"
} else {
    split_bed_dir  = "${params.split_bed}/genome"
}

// prepare split bed files
Channel
    .fromPath("${split_bed_dir}/*bed")
    .ifEmpty { exit 1, "Cannot find any bed file in directory: ${split_bed_dir}\n!" }
    .into { haplotype_beds; combine_gvcf_beds }


// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n!" }
    .into { fastqc_fq_files;  bwa_fq_files }

/*
 * FastQC
 */
process fastqc {
    tag "${name}"
    
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
    tag "${sample_name}"

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

    input:
    set sample_name, file(reads) from bwa_fq_files
    file bwa_index_file from bwa_index_file
    file genome_fa from genome_fa

    output: 
    file "${sample_name}.bam" into unsort_bam

    cpus = 16

    script:    
    """
    bwa mem -M -a \\
	    -R \"@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPI:350\\tPL:Illumina\\tCN:TCuni\" \\
	    -t ${task.cpus} \\
	    -K 10000000 \\
	    ${bwa_index_file}/${genome_fa.getName()} \\
	    ${reads[0]} \\
	    ${reads[1]} \\
	    | \\
	samtools view -O bam \\
	    --threads ${task.cpus} \\
	    -o ${sample_name}.bam
    """
}


/*
* Sort bam
*/
process sort_bam {
    tag "${sample_name}"

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

    input:
    file bam from unsort_bam

    output: 
    file "${sample_name}.sort.bam" into samtools_stats_bam, to_rmdup_bam

    cpus = 8

    script:
    sample_name = bam.baseName
    """	
    samtools fixmate --threads ${task.cpus} \\
        -m ${bam} ${sample_name}.fixmate.bam

    samtools sort -m 2400M --threads ${task.cpus} \\
	    -o ${sample_name}.sort.bam \\
	    ${sample_name}.fixmate.bam
    """
}

/*
* Reads coverage stats
*/
process reads_cov_stats {

    tag "${sample_name}"

    module "samtools/1.9"

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

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
    samtools stats \\
	    --threads ${task.cpus} \\
	    --target-regions  ${cds_bed_file} \\
	    ${samtools_stats_bam} \\
	    > ${sample_name}.cds.stat

    samtools stats \\
	    --threads ${task.cpus} \\
	    --target-regions  ${exon_bed_file} \\
	    ${samtools_stats_bam} \\
	    > ${sample_name}.exon.stat

    samtools stats \\
	    --threads ${task.cpus} \\
	    ${samtools_stats_bam} \\
	    > ${sample_name}.genome.stat 
    """
}


/*
* Reads mark duplication and recalibration
*/
process bam_remove_duplicate {
    tag "${sample_name}"

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

    input:
    file bam from to_rmdup_bam
  
    output:
    file "${sample_name}.rmdup.bam" into br_rmdup_bam
  
    cpus = 8

    script:
    sample_name = bam.baseName - '.sort'
    """
    samtools markdup -r -s \\
	    --threads ${task.cpus} \\
	    ${bam} \\
	    ${sample_name}.rmdup.bam
    """
}


/*
* bam BaseRecalibrator
*/
process bam_BaseRecalibrator {
    tag "${sample_name}"

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

    when:
    params.known_vcf    

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
    gatk BaseRecalibrator \\
        --reference ${refer} \\
        --input ${sample_name}.rmdup.bam \\
        --output ${sample_name}.recal.table \\
        --known-sites ${known_vcf}
    """
}

/*
* bam ApplyBQSR
*/
process bam_ApplyBQSR {
    tag "${sample_name}"

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

    when:
    params.known_vcf

    input:
    file bam from bqsr_rmdup_bam
    file recal_table_file from recal_table
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict
  
    output:
    file "${sample_name}.bqsr.bam" into bqsr_bam, sample_bam
  
    cpus = 8

    script:
    sample_name = bam.baseName - '.rmdup'
    """       
    gatk ApplyBQSR \\
        --bqsr-recal-file ${recal_table_file} \\
        --input ${bam} \\
        --output ${sample_name}.bqsr.bam \\
        --read-filter AmbiguousBaseReadFilter \\
        --read-filter MappingQualityReadFilter \\
        --read-filter NonZeroReferenceLengthAlignmentReadFilter \\
        --read-filter ProperlyPairedReadFilter \\
        --minimum-mapping-quality 30 
    """
}

/*
*  GATK HaplotypeCaller
*/
//saveAs: {filename -> filename.indexOf("hc.g.vcf.gz.tbi") > 0 ? null : "$filename"}
process gatk_HaplotypeCaller {
    tag "${sample_name}|${chr_name}"

    publishDir "${params.outdir}/gvcf/each_sample/${sample_name}" , mode: 'copy'

    when:
    params.known_vcf        

    input:
    file bam from bqsr_bam
    each file(bed) from haplotype_beds
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    

    output:
    file "${sample_name}.${chr_name}.hc.g.vcf.gz" into sample_gvcf
    file "${sample_name}.${chr_name}.hc.g.vcf.gz.tbi" into sample_gvcf_index
    
    cpus = 8

    script:
    sample_name = bam.baseName - '.bqsr'
    chr_name = bed.baseName
    """
    gatk HaplotypeCaller  \\
        --input ${bam} \\
        --output ${sample_name}.${chr_name}.hc.g.vcf.gz \\
        --reference ${refer} \\
        --intervals ${bed} \\
        --emit-ref-confidence GVCF \\
        --read-filter AmbiguousBaseReadFilter \\
        --read-filter MappingQualityReadFilter \\
        --read-filter NonZeroReferenceLengthAlignmentReadFilter \\
        --read-filter ProperlyPairedReadFilter \\
        --minimum-mapping-quality 30         
    """
}

/*
* GATK CombineGVCFs
*/
process gatk_CombineGVCFs {
    tag "Chrom: ${chr_name}"

    publishDir "${params.outdir}/gvcf/all_sample", mode: 'copy'

    when:
    params.known_vcf    

    input:
    file ('gvcf/*') from sample_gvcf.collect()
    file ('gvcf/*') from sample_gvcf_index.collect()
    each file(bed) from combine_gvcf_beds
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    
    
    output:
    file "all_sample.${chr_name}.g.vcf.gz" into merged_sample_gvcf
    file "all_sample.${chr_name}.g.vcf.gz.tbi" into merged_sample_gvcf_index
    
    script:
    chr_name = bed.baseName
    """
    ls gvcf/*.${chr_name}.hc.g.vcf.gz > ${chr_name}.gvcf.list

    gatk CombineGVCFs \\
    	--output all_sample.${chr_name}.g.vcf.gz \\
	    --reference ${refer} \\
	    --variant ${chr_name}.gvcf.list
    """
}

/*
* GATK GenotypeGVCFs
*/
process gatk_GenotypeGVCFs {
    tag "Chrom: ${chr_name}"

    publishDir "${params.outdir}/vcf/all_sample/", mode: 'copy'

    when:
    params.known_vcf    

    input:
    file gvcf from merged_sample_gvcf
    file gvcf_index from merged_sample_gvcf_index
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    
    
    output:
    file "${vcf_prefix}.vcf.gz" into merged_sample_vcf
    file "${vcf_prefix}.vcf.gz.tbi" into merged_sample_vcf_index
    
    script:
    vcf_prefix = gvcf.baseName - '.g.vcf'
    chr_name = vcf_prefix - 'all_sample.'
    """
    gatk GenotypeGVCFs \\
        --reference ${refer} \\
        --variant ${gvcf} \\
        --output ${vcf_prefix}.vcf.gz
    """
}

/*
* Concat vcf
*/
process concat_vcf {

    publishDir "${params.outdir}/vcf/all_chr", mode: 'copy'

    when:
    params.known_vcf    

    input:
    file ('vcf/*') from merged_sample_vcf.collect()
    file ('vcf/*') from merged_sample_vcf_index.collect()

    output:
    file "all_sample.raw.vcf.gz" into all_sample_raw_vcf
    file "all_sample.raw.vcf.gz.tbi" into all_sample_raw_vcf_idx
    
    script:
    """
    bcftools concat \\
	    vcf/*.vcf.gz | \\
        bgzip > all_sample.raw.vcf.gz

    tabix -p vcf all_sample.raw.vcf.gz
    """
}

/*
* Basic quality filter
*/
process vcf_base_qual_filter {

    publishDir "${params.outdir}/vcf/all_chr", mode: 'copy'

    when:
    params.known_vcf    

    input:
    file raw_vcf from all_sample_raw_vcf
    file raw_vcf_idx from all_sample_raw_vcf_idx
    
    output:
    file "all_sample.hq.vcf.gz" into all_hq_vcf, all_hq_vcf_for_extract
    file "all_sample.hq.vcf.gz.tbi" into all_hq_vcf_idx, all_hq_vcf_idx_for_extract
    
    script:
    """
    bcftools filter -s LowQual -e '%QUAL<${params.quality} || INFO/DP<${params.depth}' \\
	    all_sample.raw.vcf.gz | \\
        grep -v LowQual | bgzip > all_sample.hq.vcf.gz

    tabix -p vcf all_sample.hq.vcf.gz
    """
}

/*
* snpeff for combined vcf
*/
process snpEff_for_all {

    publishDir "${params.outdir}/vcf/all_chr", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db

    input:
    file vcf from all_hq_vcf
    file vcf_idx from all_hq_vcf_idx
    
    output:
    file "all_sample.hq.vcf.stat.csv"
    file "all_sample.hq.vcf.stat.html" 
    file "all_sample.ann.vcf.gz" into all_sample_anno_vcf
    file "all_sample.ann.vcf.gz.tbi" into all_sample_anno_vcf_idx
    
    script:
    """
    java -Xmx10g -jar ${params.snpEff}/snpEff.jar \\
        -c ${params.snpEff}/snpEff.config \\
        -csvStats all_sample.hq.vcf.stat.csv \\
        -htmlStats all_sample.hq.vcf.stat.html \\
        -v ${params.snpEff_db}  \\
        ${vcf} \\
        | bgzip > all_sample.ann.vcf.gz
    
    tabix -p vcf all_sample.ann.vcf.gz
    """
}


/*
* extract vcf file for each sample
*/
process extract_sample_vcf {
    tag "${sample_name}"

    publishDir "${params.outdir}/vcf/each_sample/${sample_name}", mode: 'copy'

    when:
    params.known_vcf    

    input:
    file vcf from all_hq_vcf_for_extract
    file vcf_idx from all_hq_vcf_idx_for_extract
    file bam from sample_bam

    output:
    file "${sample_name}.hq.vcf.gz" into sample_hq_vcf
    file "${sample_name}.hq.vcf.gz.tbi" into sample_hq_vcf_idx
    
    script:
    sample_name = bam.baseName - '.bqsr'
    """
    bcftools view -Ov -s ${sample_name} ${vcf} | \\
	    bcftools filter -i 'GT="alt"' | \\
        bgzip > ${sample_name}.hq.vcf.gz

    tabix -p vcf ${sample_name}.hq.vcf.gz
    """
}

/*
* snpEff for each sample
*/
process snpEff_for_sample {

    tag "${sample_name}"

    publishDir "${params.outdir}/vcf/each_sample/${sample_name}", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db

    input:
    file vcf from sample_hq_vcf
    file vcf_idx from sample_hq_vcf_idx
    
    output:
    file "${sample_name}.hq.vcf.stat.csv"
    file "${sample_name}.hq.vcf.stat.html"
    file "${sample_name}.ann.vcf.gz" into anno_vcf
    file "${sample_name}.ann.vcf.gz.tbi"
    
    script:
    sample_name = vcf.baseName - '.hq.vcf'
    """
    java -Xmx10g -jar ${params.snpEff}/snpEff.jar \\
        -c ${params.snpEff}/snpEff.config \\
        -csvStats ${sample_name}.hq.vcf.stat.csv \\
        -htmlStats ${sample_name}.hq.vcf.stat.html \\
        -v ${params.snpEff_db}  \\
        ${vcf} \\
        | bgzip > ${sample_name}.ann.vcf.gz

    tabix -p vcf ${sample_name}.ann.vcf.gz
    """
}

/*
* SNP Table
*/
process snp_table {

    publishDir "${params.outdir}/vcf/all_chr", mode: 'copy'

    input:
    file vcf from all_sample_anno_vcf
    file vcf_idx from all_sample_anno_vcf_idx
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict 

    output:
    file "all_sample_gatk.table.txt" into gatk_vcf_table
    file "all_sample.vcf.table.txt" into om_vcf_table

    script:
    """
    java -jar /public/software/GATK/GATK.3.8/GenomeAnalysisTK.jar \\
        -T VariantsToTable \\
        -R ${refer} \\
        -V ${vcf} \\
        -F CHROM -F POS -F REF -F ALT \\
        -GF AD -GF DP -GF GQ -GF PL \\
        -o all_sample_gatk.table.txt
    
    gunzip -c ${vcf} > ${vcf.baseName}

    python /public/scripts/Reseq/omtools/extractTableFromsnpEff.py \\
        -v ${vcf.baseName} \\
        -o all_sample.vcf.table.txt
    """
    }