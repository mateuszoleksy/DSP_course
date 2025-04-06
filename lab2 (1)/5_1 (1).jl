##
function quantize(x, l)
    min = round(abs(l[1]-x))
    indexid = 1
    for i in 1:1:length(l)
        if (min > round(abs(l[i]-x)))
            min = round(abs(l[i]-x))
            indexid = i
        end
    end
    return l[indexid]
end


x = range(1,10,100)
b = sin.(x).*4
y = Float64[]
push!(y, 0)
push!(y, 1)
push!(y, 2)
push!(y, 3)
push!(y, 4)
push!(y, -1)
push!(y, -2)
push!(y, -3)
push!(y, -4)
using CairoMakie
matrix = Float64[]
for i in 1:1:length(x)
push!(matrix, quantize(b[i], y))
end
f = lines(x, matrix)
f

##