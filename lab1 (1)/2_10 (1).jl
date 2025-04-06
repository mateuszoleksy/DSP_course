##
function sawtooth_wave(offset, time)
    x=Float64[]
    y=Float64[]
    push!(x,offset)
    push!(y,0)
    for i in 1:1:time*100
            push!(x, i/100+offset)
        end
    for i in 1:1:time
        for k in 1:1:100
            push!(y, k*(-1)/100+0.5)
        end
    end
    return [x,y]
end

using CairoMakie
f = Figure()
matrix = sawtooth_wave(2, 10)
f = lines(matrix[1], matrix[2])
f
##