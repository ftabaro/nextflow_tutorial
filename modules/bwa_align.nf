/*
========================================================================================
   BWA-ALIGN module
========================================================================================
   Website: http://bio-bwa.sourceforge.net/
========================================================================================
*/

// Parameter definitions
params.CONTAINER = "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:66ed1b38d280722529bb8a0167b0cf02f8a0b488-0"
params.OUTPUT = "bwa_align"

process BWA_ALIGN {
 
 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )
   
 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )
 
 // show in the log which input file is analysed
 tag( "${sample_id}" )
 
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