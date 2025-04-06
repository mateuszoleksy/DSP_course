
##
function signal(n)
    return 2*sin(2*0.04*n + pi/4)
end

using CairoMakie
f = Figure()

wektor=Float64[]
for i in 1:1:256
push!(wektor, signal(i/1000))
end
lines(range(1,256,256), wektor)
f
@show wektor
##  