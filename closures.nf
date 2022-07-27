square = { it * it }
x = [ 1, 2, 3, 4 ]
println(x)
println(x.collect(square))

prefix = {"chr${it}"}
x = x.collect(prefix)
println x
