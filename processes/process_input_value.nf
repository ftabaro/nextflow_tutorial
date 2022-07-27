//process_input_file.nf
nextflow.enable.dsl=2

/*
 * Index the reference genome for use by bwa and samtools.
 */
process BWA_INDEX {

  input:
  path genome

  script:
  """
  bwa index ${genome}
  """
}

ref_ch = Channel.fromPath("/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta")  

workflow {
  BWA_INDEX( ref_ch )
}
