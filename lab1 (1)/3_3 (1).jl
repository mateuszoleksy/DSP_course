function energy(x)
    return sum(abs.(x)).^2
end

x = range(1,100,100)
y = sin.(x)

@show energy(y)