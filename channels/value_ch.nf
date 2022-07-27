//Creates an empty value channel.
ch1 = Channel.value( 'GRCh38' )
//Creates a value channel and binds a string to it.
ch2 = Channel.value( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'] )
//Creates a value channel and binds a list object to it that will be emitted as a single item.
ch3 = Channel.value( [chr1 : "248956422", 'chr2' : 242193529, 'chr3' : 198295559] )
//The value method can only take 1 argument, however, this can be a single list containing several elements.
ch1.view()
ch2.view()
ch3.view()

