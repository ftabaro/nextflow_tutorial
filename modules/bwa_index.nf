/*
========================================================================================
   BWA-INDEX module
========================================================================================
   Website: http://www.htslib.org/doc/samtools.html
========================================================================================
*/

// Parameter definitions
params.CONTAINER = "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"
params.OUTPUT = "bwa_index"

process BWA_INDEX {
 
 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )
   
 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )
 
 // show in the log which input file is analysed
 tag( "${genome}" )
 
 input:
 path( genome )

 output:
 tuple path( genome ), path( "*" ), emit: bwa_index

 script:
 """
 bwa index ${genome} 
 """
}