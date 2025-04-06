##
function rms(x)
    return sqrt(sum(abs.(x))^2/length(x))
end

x = range(1,100,100)
y = sin.(x)

@show rms(y)
##