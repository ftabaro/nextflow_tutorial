read_pair_ch = Channel.fromFilePairs("/workspace/nextflow_tutorial/data/untrimmed_fastq/*_{1,2}.fastq.gz", checkIfExists: true)
read_pair_ch.view()
