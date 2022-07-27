aligner_list = ['salmon', 'kallisto']

aligner_ch = Channel.fromList(aligner_list)

aligner_ch.view()
