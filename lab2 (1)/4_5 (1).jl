##
function interpolate(t, m, s, kernel::Function=sinc)
    y = Float64[]

    for i in 1:1:length(m)-1
        push!(y, s[i] + kernel((s[i+1]-s[i])/t*(m[i+1]-m[i])))
    end
    push!(y, y[end])
    
    return y
end


x = range(1,20,100)
y = sin.(x)
using CairoMakie
f = lines(x, y)
f
y = interpolate(1, x, y)
f = lines(x, y)
f

##