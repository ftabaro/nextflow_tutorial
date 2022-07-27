//process_publishDir_semantic.nf
nextflow.enable.dsl=2

process FASTQC {

    tag {"FASTQC $sample_id"}
    label 'process_low'
    cpus 2

    publishDir "results/fastqc_html", pattern: "*.html", mode: 'copy'
    publishDir "results/fastqc_zip", pattern: "*.zip", mode: 'copy'
    
    input:
    tuple val( sample_id ), path( reads )

    output:
    tuple val( sample_id ), path( "*.html" ), path( "*.zip" )

    script:
    """
    fastqc ${reads}
    """

}

reads_ch = Channel.fromFilePairs("/workspace/nextflow_tutorial/data/trimmed_fastq/SRR2584863_{1,2}.trim.fastq.gz", checkIfExists: true)

workflow {
  
  FASTQC( reads_ch )
  
}
