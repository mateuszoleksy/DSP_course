function power(x)
    return (sum(abs.(x)).^2)./length(x)
end

x = range(1,100,100)
y = sin.(x)

@show power(y)