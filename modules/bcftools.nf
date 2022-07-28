/*
========================================================================================
   BCFTOOLS module
========================================================================================
   Website: http://www.htslib.org/doc/bcftools.html
========================================================================================
*/

// Parameter definitions
params.CONTAINER = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
params.OUTPUT = "vcf"

process BCFTOOLS_MPILEUP {

 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )

 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )

 // show in the log which input file is analysed
 tag( "${sample_id}" )
 
 input:
 tuple path( genome ), path( "*" ), val( sample_id ), path( sortedbam ), path("*")

 output:
 tuple val( sample_id ), path( "${sample_id}_raw.bcf" ), emit: raw_bcf

 script:
 """
 bcftools mpileup -O b -o "${sample_id}_raw.bcf" -f ${genome} ${sortedbam}
 """
}

/*
* Detect the single nucleotide variants (SNVs).
*/
process BCFTOOLS_CALL {
 
 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )

 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )

 // show in the log which input file is analysed
 tag( "${sample_id}" )

 input:
 tuple val( sample_id ), path( rawbcf )

 output:
 tuple val( sample_id ), path( "${sample_id}_variants.vcf" ), emit: variants_vcf

 script:
 """
 bcftools call --ploidy 1 -m -v -o "${sample_id}_variants.vcf" ${rawbcf}
 """
}

process VCFUTILS {
 // where to store the results and in which way
 publishDir( params.OUTPUT, mode: 'copy' )

 // indicates to use as container the value indicated in the parameter
 container( params.CONTAINER )

 // show in the log which input file is analysed
 tag( "${sample_id}" )

 input:
 tuple val( sample_id ), path( rawvcf )

 output:
 tuple val( sample_id ), path( "${sample_id}_final_variants.vcf" ), emit: final_variants_vcf

 script:
 """
 vcfutils.pl varFilter ${rawvcf} > "${sample_id}_final_variants.vcf"
 """
}