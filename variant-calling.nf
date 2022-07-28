
/*
========================================================================================
    Variant-Calling Nextflow Workflow
========================================================================================
    Github   : 
    Contact  :     
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

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
    Include Modules
========================================================================================
*/

include { FASTQC }                                    from "./modules/fastqc" addParams(OUTPUT: "${params.outdir}/fastqc")
include { BWA_INDEX  }                                from "./modules/bwa_index" addParams(OUTPUT: "${params.outdir}/bwa_index")
include { BWA_ALIGN  }                                from "./modules/bwa_align" addParams(OUTPUT: "${params.outdir}/bwa_align")
include { SAMTOOLS_SORT; SAMTOOLS_INDEX }             from "./modules/samtools" addParams(OUTPUT: "${params.outdir}/sorted_bam")
include { BCFTOOLS_MPILEUP; BCFTOOLS_CALL; VCFUTILS } from "./modules/bcftools" addParams(OUTPUT: "${params.outdir}/vcf")

/*
========================================================================================
    Create Channels
========================================================================================
*/

ref_ch = Channel.fromPath( params.genome, checkIfExists: true  )  
reads_ch = Channel.fromFilePairs( params.reads, checkIfExists: true ) 

/*
========================================================================================
    WORKFLOW - Variant Calling
========================================================================================
*/

workflow {

    FASTQC( reads_ch )
    BWA_INDEX( ref_ch )
    BWA_ALIGN( BWA_INDEX.out.bwa_index.combine(reads_ch) )
    SAMTOOLS_SORT( BWA_ALIGN.out.aligned_bam )
    SAMTOOLS_INDEX( SAMTOOLS_SORT.out.sorted_bam )
    BCFTOOLS_MPILEUP( BWA_INDEX.out.bwa_index.combine(SAMTOOLS_INDEX.out.aligned_sorted_bam) )
    BCFTOOLS_CALL( BCFTOOLS_MPILEUP.out.raw_bcf )
    VCFUTILS( BCFTOOLS_CALL.out.variants_vcf )

}

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

