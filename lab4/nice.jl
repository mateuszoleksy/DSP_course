function fft3(x::Vector{ComplexF64})::Vector{ComplexF64}
    matrix = ComplexF64[]
    suma1 = 0
    suma2 = 0  
    N = length(x)
    len = N/2-1
    
    for k in 1:1:Int(len)
      angle2 = -2*pi*k/N
          for h in 1:1:Int(len)
              angle = -2*pi*h*k/(N/2)
              suma1 += x[h*2]*exp(1im * angle)
              suma2 += x[h*2+1]*exp(1im * angle)
        end
      push!(matrix, suma1 + exp(1im * angle2).*suma2)
      suma1 = 0
      suma2 = 0
 end
  return matrix
end

N = 1
fp = 100
test = sin.(40*pi*range(0,N,fp)) +  sin.(10*pi*range(0,N,fp)) .+ 0im
y = fft3(test)