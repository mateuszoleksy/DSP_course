function mean(x)
    return sum(abs.(x))/length(x)
end

x = range(1,100,100)
y = sin.(x)

@show mean(y)