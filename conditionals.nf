x = 12
if( x > 10 )
    println "$x is greater the 10"

// Null, empty strings and empty collections are evaluated to false in groovy

list = []

if( list )
  println list
else
  println 'The list is empty'
