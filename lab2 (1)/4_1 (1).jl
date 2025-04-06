
##
function signal(n)
    return 2*sin(200*n)
end

using CairoMakie
f = Figure()

wektor=Float64[]
for i in 1:1:1000
push!(wektor, signal(i*4/1000))
end
lines(range(1,1000,1000), wektor)
f
@show wektor
##  