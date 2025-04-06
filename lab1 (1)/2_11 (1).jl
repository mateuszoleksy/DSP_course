##
function triangular_wave(offset, time)
    x=Float64[]
    y=Float64[]
    push!(x,offset)
    push!(y,0)
    for i in 1:1:time*400+100
            push!(x, i/100+offset)
        end
    for k in 1:1:100
        push!(y, k/100)
    end
    for i in 1:1:time
        for k in 1:1:200
            push!(y, k/100*(-1)+1)
        end
        for k in 1:1:200
            push!(y, k/100-1)
        end
    end
    return [x,y]
end

using CairoMakie
f = Figure()
matrix = triangular_wave(2, 10)
lines(matrix[1], matrix[2])
##