function impuse_repeater(g::Function, t1::Real, t2::Real)

matrix = Float64[]
append!(matrix,  g.(range(t1,t2,100)))
append!(matrix, g.(range(t1,t2,100)))
append!(matrix, g.(range(t1,t2,100)))
append!(matrix, g.(range(t1,t2,100)))

end

x = impuse_repeater(cos, 4, 2*3)
using CairoMakie
f = Figure()
f = lines(range(3,24,400), x)
f