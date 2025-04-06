##



function fft2(x)
    matrix = ComplexF64[]
    suma1 = 0
    suma2 = 0  
    N = length(x)
    len::Integer = N/2-1
    
    for k in 1:1:N::Integer
      angle2 = -2*pi*k/N
          for h in 1:1:len::Integer
              angle = -2*pi*h*k/(N/2)
              suma1 += x[h*2]*exp(1im * angle)
              suma2 += x[h*2+1]*exp(1im * angle)
        end
      push!(matrix, suma1 + exp(1im * angle2).*suma2)
      suma1 = 0
      suma2 = 0
 end
  return abs.(matrix)
end

N = 1
fp = 100
y = sin.(40*pi*range(0,N,fp)) +  sin.(10*pi*range(0,N,fp))
T = 1/fp
freq = range(0,N,fp)/T
using CairoMakie
f = lines(freq, fft2(y))
##