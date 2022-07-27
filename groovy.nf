println("Hello, World!")

// Hello! 

/*
Block
*/

//int − This is used to represent whole numbers.
my_int = 1

//float − This is used to represent  floating point numbers.
my_float = 3.1499392

//Boolean − This represents a Boolean value which can either be true or false.
my_bool = false

//String - These are text literals which are represented in the form of chain of characters
my_string = "chr1"

// A block of text that span multiple lines can be defined by delimiting it with triple single `'''` or double quotes `"""`:
my_text = """
          This is a multi-line
          using triple quotes.
          """

// To display the value of a variable to the screen in Groovy, we can use the `println` method passing the variable name are a parameter.
println(my_int)
println(my_float)
println(my_bool)
println(my_string)
println(my_text)

// String Interpolation
// To use a variable inside a single or multi-line double quoted string "" prefix the variable name with a $ to show it should be interpolated. 
println("processing chromosome $my_int")
println("value of pi is $my_float")

def x = 'local_variable_def'
println(x)

kmers = [11,21,27,31]
// You can access a given item in the list with square-bracket notation []. These positions are numbered starting at 0, so the first element has an index of 0.
println(kmers[0])
// Lists can also be indexed with negative indexes
println(kmers[-1])
// The first three elements Lists elements using a range.
println(kmers[0..2])
// String interpolation for Lists - To use an expression like `kmer[0..2]` inside a double quoted String `""` we use the `${expression}` syntax, similar to Bash/shell scripts
println("The first three elements in the Lists are: ${kmers[0..2]}")

// To get the length of the list
println(kmers.size())

//inside a string need we need to use the ${} syntax
println("list size is:  ${kmers.size()}")

// To retrieve items in a list
println kmers.get(1)




