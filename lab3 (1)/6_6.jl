##
#program do obliczania dft
function rdft(x::AbstractVector)::Vector
    matrix = Float64[]
    suma = 0
    for k in 1:1:length(x)
        for i in 1:1:length(x)
           angle = 2*pi*i*k/length(x)
           suma += x[i]*exp(-1im * angle)
           #lub k*i%length arrows 
        end
        push!(matrix, real(suma))
        suma = 0
    end
    return matrix
end

x = [20, 5]
y = sin.(x)

@show rdft(y)
using CairoMakie
f = scatter(x,abs.(rdft(y)))
f
##