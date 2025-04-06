##
function rfft(N::Integer, fs::Integer)
    x = sin.(2*pi*range(0,N,N))+cos.(2*2*pi*range(0,N,N))
    matrix = ComplexF64[]
    suma = 0
    suma2 = 0
    len::Integer = N/2-1
        for k in 1:1:len::Integer
             for h in 1:1:len::Integer
               angle = -2*pi*h*k/len
               angle2 = -2*pi*k/len
               suma += x[h]*exp(-1im * angle)
               suma2 += exp(-1im * angle2)
               #lub k*i%length arrows 
            end
            push!(matrix, suma*(1im) + suma*1im/(suma2*1im))
            suma = 0
            suma2 = 0
     end
     len = N/2
     for k in len::Integer:1:N::Integer
            for h in len::Integer:1:N-1
              angle = -2*pi*h*(k-N/2)/len
              angle2 = -2*pi*k/len
              suma += x[h]*exp(-1im * angle)
              suma2 += exp(-1im * angle2)
              #lub k*i%length arrows 
          end
          push!(matrix, suma*(1im) + suma*1im/(suma2*1im))
          suma = 0
          suma2 = 0
 end
  return abs.(matrix*1im/N)
end

@show rfft(1000, 100)
using CairoMakie
x = plot(range(0,1000,1000), rfft(1000, 100))
##