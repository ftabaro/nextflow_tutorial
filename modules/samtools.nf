/*
========================================================================================
   SAMTOOLS module
========================================================================================
   Website: http://www.htslib.org/doc/samtools.html
========================================================================================
*/

// Parameter definitions
params.CONTAINER = "quay.io/biocontainers/samtools:1.14--hb421002_0"
params.OUTPUT = "sorted_bam"

process SAMTOOLS_SORT {

 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )

 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )

 // show in the log which input file is analysed
 tag( "${sample_id}" )

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
 
 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )

 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )

 // show in the log which input file is analysed
 tag( "${sample_id}" )

 input:
 tuple val( sample_id ), path( sortedbam )

 output:
 tuple val( sample_id ), path( "${sample_id}.aligned.sorted.bam" ), path("*"), emit: aligned_sorted_bam

 script:
 """
 samtools index ${sortedbam}
 """
}