/*
========================================================================================
   Variant-Calling Nextflow Workflow
========================================================================================
   Github   : 
   Contact  :     
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Pipeline Input parameters

params.outdir = 'results'
params.genome = "/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta"
params.reads = "/workspace/nextflow_tutorial/data/trimmed_fastq/*_{1,2}.trim.fastq.gz"

println """\
         V A R I A N T-C A L L I N G - N F   P I P E L I N E
         ===================================
         genome       : ${params.genome}
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

/*
========================================================================================
   Create Channels
========================================================================================
*/

ref_ch = Channel.fromPath( params.genome, checkIfExists: true )  
reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 

/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

    FASTQC( reads_ch )
    BWA_INDEX( ref_ch )
    BWA_ALIGN( BWA_INDEX.out.bwa_index.combine(reads_ch) ) // https://www.nextflow.io/docs/latest/process.html#understand-how-multiple-input-channels-work
    SAMTOOLS_SORT( BWA_ALIGN.out.aligned_bam )
    // Enter the rest of the processes for variant calling based on the bash script below
    SAMTOOLS_INDEX( SAMTOOLS_SORT.out.sorted_bam )
    BCFTOOLS_MPILEUP( BWA_INDEX.out.bwa_index.combine(SAMTOOLS_SORT.out.sorted_bam) )
    BCFTOOLS_CALL( BCFTOOLS_MPILEUP.out.raw_bcf )
    VCFUTILS( BCFTOOLS_CALL.out.variants )
}

/*
========================================================================================
   Processes
========================================================================================
*/

/*
 * Align reads to reference genome & create BAM file.
 */
process FASTQC {
    tag{"FASTQC ${reads}"}
    label 'process_low'

    publishDir("${params.outdir}/fastqc_trim", mode: 'copy')
    
    input:
    tuple val( sample_id ), path( reads )

    output:
    path( "*_fastqc*" ), emit: fastqc_out

    script:
    """
    fastqc ${reads}
    """
}

/*
 * Index the reference genome for use by bwa and samtools.
 */
process BWA_INDEX {
  tag{"BWA_INDEX ${genome}"}
  label 'process_low'

  publishDir("${params.outdir}/bwa_index", mode: 'copy')
  
  input:
  path genome

  output:
  tuple path( genome ), path( "*" ), emit: bwa_index

  script:
  """
  bwa index ${genome} 
  """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BWA_ALIGN {
    tag{"BWA_ALIGN ${sample_id}"}
    label 'process_medium'

    publishDir("${params.outdir}/bwa_align", mode: 'copy')
    
    input:
    tuple path( genome ), path( "*" ), val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.bam" ), emit: aligned_bam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem \$INDEX ${reads} > ${sample_id}.aligned.sam
    samtools view -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
    """
}

/*
 * Convert the format of the alignment to sorted BAM.
 */
process SAMTOOLS_SORT {
  tag{"SAMTOOLS_SORT ${sample_id}"}
  label 'process_low'

  publishDir("${params.outdir}/bam_align", mode: 'copy')

  input:
  tuple val( sample_id ), path( bam )

  output:
  tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), emit: sorted_bam

  script:
  """
  samtools sort -o "${sample_id}.aligned.sorted.bam" ${bam}
  """
}

/*
 * Index the BAM file for visualization purpose
 */
process SAMTOOLS_INDEX {
  tag{"SAMTOOLS_INDEX ${sample_id}"}
  label 'process_low'

  publishDir("${params.outdir}/bam_align", mode: 'copy')

  input:
  tuple val( sample_id ), path( bam )

  output:
  tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam.bai" ), emit: sorted_bam

  script:
  """
  samtools index ${sample_id}.aligned.sorted.bam
  """
}

/*
 * Calculate the read coverage of positions in the genome.
 */
process BCFTOOLS_MPILEUP {
  tag{"BCFTOOLS_MPILEUP ${sample_id}"}
  label 'process_low'

  publishDir("${params.outdir}/bcf", mode: 'copy')

  input:
  tuple path( genome ), path( "*" ), val( sample_id ), path( sorted_bam )

  output:
  tuple val( sample_id ), path( "${sample_id}_raw.bcf" ), emit: raw_bcf

  script:
  """
  bcftools mpileup -O b -o ${sample_id}_raw.bcf -f ${genome} ${sorted_bam}
  """
}

/*
 * Detect the single nucleotide variants (SNVs).
 */
process BCFTOOLS_CALL {
  tag{"BCFTOOLS_CALL ${sample_id}"}
  label 'process_low'

  publishDir("${params.outdir}/variant_calls", mode: 'copy')

  input:
  tuple val(sample_id), path( raw_bcf )

  output:
  tuple val( sample_id ), path( "${sample_id}_variants.vcf" ), emit: variants

  script:
  """
  bcftools call --ploidy 1 -m -v -o ${sample_id}_variants.vcf ${raw_bcf}
  """
}

/*
 * Filter and report the SNVs in VCF (variant calling format).
 */
process VCFUTILS {
  tag{"BCFTOOLS_CALL ${sample_id}"}
  label 'process_low'

  publishDir("${params.outdir}/vcf", mode: 'copy')

  input:
  tuple val( sample_id ), path( variants )

  output:
  tuple val( sample_id ), path( "${sample_id}_final_variants.vcf" ), emit: final_variants

  script:
  """
  vcfutils.pl varFilter ${variants} > ${sample_id}_final_variants.vcf
  """
}

/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/
