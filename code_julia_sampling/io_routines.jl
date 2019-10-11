
function writedata(data::Array, filename::AbstractString)
	iostream = open(filename, "w")
	writedlm(iostream, data)
	close(iostream)
end

function readvec(file::AbstractString)
	iostream = open(file, "r")
	vec1 = readdlm(iostream, comments=true, comment_char='#')[:,1]
	close(iostream)
	return vec1
end

# Significant digits
function sig(var, sigdig)
	return round(var,sigdigits=sigdig)
end