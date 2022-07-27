#!/usr/bin/env nextflow
nextflow.enable.dsl=2

chromosome_ch = Channel.of( 'chr1','chr3','chr5','chr7' )
chromosome_ch.bind(1..22, 'X', 'Y')
chromosome_ch.view()
