//process_tuple_input.nf
nextflow.enable.dsl=2

process TUPLEINPUT{
  input:
  tuple val(sample_id), path(reads)
  
  script:
  """
  echo ${sample_id}
  echo ${reads}
  """
}

reads_ch = Channel.fromFilePairs("/workspace/nextflow_tutorial/data/trimmed_fastq/SRR2584863_{1,2}.trim.fastq.gz")

workflow {
  
  TUPLEINPUT(reads_ch)
  
}
