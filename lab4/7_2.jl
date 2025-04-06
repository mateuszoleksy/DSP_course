##
function fft(x)
    matrix = ComplexF64[]
    suma1 = 0
    suma2 = 0
    N = length(y)
    len::Integer = N/2-1
    
    for k in 1:1:len::Integer
      angle2 = -2*pi*k/len
          for h in 1:1:len::Integer
              angle = -2*pi*h*k/len
               
              suma1 += x[h]*exp(1im * angle)
      end
      push!(matrix, suma1 + exp(1im * angle2)*suma1)
 end
 suma1 = 0
 len = N/2
 for k in len:1:N::Integer
    angle2 = -2*pi*k/len
        for h in len:1:N::Integer
            angle = -2*pi*h*k/len
            
            suma1 += x[h]*exp(1im * angle)
    end
    push!(matrix, suma1 + exp(1im * angle2)*suma1)
end
  return abs.(matrix)
end

N = 10000
fp = 100
y = sin.(2*pi*range(0,N,fp))
T = 1/fp
freq = range(0,N,fp)/T

using CairoMakie
x = scatter(freq, fft(y))
##