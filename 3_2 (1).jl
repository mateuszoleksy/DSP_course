function peak2peak(matrix)

    max = matrix[begin]
    min = matrix[begin]
    for i in 2:1:length(matrix)
        if (max < matrix[i])
            max = matrix[i]
        end
        if (min > matrix[i])
            min = matrix[i]
        end
    end
    @show min
    @show max
    return abs(max-min)

end

wektor = Float64[]
for i in 1:1:1000
push!(wektor, 4*sin(2*pi*i))
end

@show wektor

@show peak2peak(wektor)