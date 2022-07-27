//workflow.nf
nextflow.enable.dsl=2

// Initialize required parameters
params.outdir = 'results'
params.genome = "/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta"
params.reads = "/workspace/nextflow_tutorial/data/trimmed_fastq/SRR2584863_{1,2}.trim.fastq.gz"

workflow {
    
    // Create channel from path for Reference Genome
    ref_ch = Channel.fromPath( params.genome, checkIfExists: true )
    // Create channel from file-pairs for Input Fastq files
    reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 
    
    //index process takes 1 input channel as a argument
    BWA_INDEX( ref_ch )

    //bwa align process takes 2 input channels as arguments
    BWA_ALIGN( BWA_INDEX.out, reads_ch )
    
}

process BWA_INDEX {
  tag {"BWA_INDEX ${genome}"}
  label 'process_low'

  publishDir "${params.outdir}/bwa_index", mode: 'copy'
  
  input:
  path( genome )

  output:
  tuple path( genome ), path( "*" )

  script:
  """
  bwa index ${genome}
  """
}

/*
 * Align reads to reference genome & create BAM file.
 */
process BWA_ALIGN {
    tag {"BWA_ALIGN ${sample_id}"}
    label 'process_medium'

    publishDir "${params.outdir}/bwa_align", mode: 'copy'
    
    input:
    tuple path( genome ), path( "*" )
    tuple val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "${sample_id}.aligned.bam" )

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem \$INDEX ${reads} > ${sample_id}.aligned.sam
    samtools view -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
    """
}
