
##
function signal(s, i)
    return acos(s)/i
end

using CairoMakie
f = Figure()
fp = 250

@show signal(sin(200*pi*1*(1/fp)), 1*(1/fp))
##  