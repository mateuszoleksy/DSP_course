##

function noise(power, n)
    matrix=Float64[]
    for i in 1:1:n
        push!(matrix, (sqrt(power)*cos(i)-sqrt(power)*sin(i)))
    end 
return matrix
end
using CairoMakie

f = Figure()
lines(range(1,100,1000), sin.(range(1,100,1000).+noise(0.25,1000)))
f
@show noise(0.25, 1000)
##