//process_index.nf
nextflow.enable.dsl=2

process BWA_INDEX {

  tag {"BWA_INDEX $genome"}
  label 'process_low'
  cpus 1

  input:
  path genome

  output:
  path("*")

  script:
  """
  bwa index ${genome}
  """
}

ref_ch = Channel.fromPath("/workspace/nextflow_tutorial/data/ref_genome/ecoli_rel606.fasta")  

workflow {
  
  BWA_INDEX( ref_ch )
  
}
