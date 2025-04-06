
##
function peak2peak(matrix)
    max = matrix[1]
    min = matrix[1]
    for i in 2:1:length(matrix)
        if (max < matrix[i])
            max = matrix[i]
        end
        if (min > matrix[i])
            min = matrix[i]
        end
    end
    return abs(max-min)

end


wektor = Float64[]
for i in 1:1:1000
push!(wektor, abs(4*sin(2*i)+2im*cos(2*i)))
end

f = Figure()
lines(range(1,1000,1000),wektor)
@show peak2peak(wektor)
##