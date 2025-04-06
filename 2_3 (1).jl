
function noise(power, n)
    matrix=Float64[]
    for i in 1:1:n
        push!(matrix, (sqrt(power)*cos(i*pi)-sqrt(power)*sin(i*pi)))
    end 
return matrix
end

@show noise(0.25, 1000)