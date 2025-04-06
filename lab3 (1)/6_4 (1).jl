##
#program do obliczania dft
function dft(x::AbstractVector)::Vector
    matrix = []
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
    for i in length(x):-1:1
    push!(matrix, matrix[i])
    end
    return matrix
end

x = [20, 5]
y = sin.(x)
x = [-20,-5, 5, 20]

using CairoMakie
f = scatter(x,abs.(dft(y)))
f
##