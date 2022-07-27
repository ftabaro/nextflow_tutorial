roi = [ chromosome : "chr17", start: 7640755, end: 7718054, genes: ['ATP1B2','TP53','WRAP53']]

// Maps can be accessed in a conventional square-bracket syntax or as if the key was a property of the map or using the dot notation. **Note: When retrieving a value the key value is enclosed in quotes.**
println(roi['chromosome'])

//Use a dot notation            
println(roi.start)

//Use of get method                      
println(roi.get('genes'))      