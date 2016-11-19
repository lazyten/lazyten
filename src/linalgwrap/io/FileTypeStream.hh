// use operator << and format-like flags to insert data.
//
// << label(string)	start a labelled record
// 			check that no record is undetermined
// << endr		end a record (maybe this can be made implicit
// 			and is not needed explicitly)
// << double		insert a double value
// << std::vector<double> insert a list of double values
// << matrix		insert a matrix
// << vector		insert a vector
// << string		insert a comment
// << verbatim << string insert a string verbatim
//
// Printing state:
// << precision		modify precision
// << fixed		modify fixed vs scientific
// << scientific	modify fixed vs scientific
// 	collect this state inside a struct and pass
// 	that on to FileType_i objects
//
// .substream("label")  get a substream of given label (i.e. go one level lower
// 			in yaml, json, ...)
