##

using WAV

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

y, fs = wavread("H:/My Drive/CPS/cps/MOleksy (1)/lab4/fail_trombone.wav")

N = 1
T = 1/Int(fs)
freq = range(0,N,Int(fs))/T
using CairoMakie
f = lines(freq, fft2(y[1:16000]))
##