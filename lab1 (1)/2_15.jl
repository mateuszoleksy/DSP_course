##
#todo ograniczenie pasma
function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    x=Float64[]
    y=Float64[]
    push!(x,t)
    push!(y,0)
    for i in 1:1:10*100
            push!(x, i/100+t)
        end
    for i in 1:1:10
        for k in 1:1:100
            push!(y, k*(-A/T)/100+A)
        end
    end
    return [x,y]
end

using CairoMakie
f = Figure()
matrix = ramp_wave_bl(2; A=1, T=1, band=20)
f = lines(matrix[1], matrix[2])
f
##