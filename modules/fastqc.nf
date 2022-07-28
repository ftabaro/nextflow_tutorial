/*
========================================================================================
   FASTQC module
========================================================================================
   Website: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
========================================================================================
*/

// Parameter definitions
params.CONTAINER = "quay.io/biocontainers/fastqc:0.11.9--0"
params.OUTPUT = "trim_fastqc_output"

process FASTQC {
   
   // where to store the results and in which way
   publishDir( params.OUTPUT, mode: 'copy' ) 

   // indicates to use as container the value indicated in the parameter
   container( params.CONTAINER )

   // show in the log which input file is analysed
   tag( "${reads}" )
   
   input:
   tuple val( sample_id ), path( reads )

   output:
   path( "*_fastqc*" ), emit: fastqc_out

   script:
   """
   fastqc ${reads}
   """
}