##
function kronecker(n)
    x=Float64[]
    y=Float64[]
    if (n < 0)
        push!(x,2*n)
        push!(y,0)
        for i in n*2:1:0
            if (i == n)
                push!(x, i-0.0000001)
                push!(x, i)
                push!(x, i+0.0000001)
            else
                push!(x, i)
            end
        end
        for i in 2*n:1:0
            if (i == n)
                push!(y, 0)
                push!(y, 1)
                push!(y, 0)
            else
                push!(y, 0)
            end
        end
    end
    if (n > 0)
        push!(x,0)
        push!(y,0)
        for i in 0:1:2*n
            if (i == n)
                push!(x, i-0.000001)
                push!(x, i)
                push!(x, i+0.000001)
            else
                push!(x, i)
            end
        end
        for i in 0:1:2*n
            if (i == n)
                push!(y, 0)
                push!(y, 1)
                push!(y, 0)
            else
                push!(y, 0)
            end
        end
    end
    if (n == 0)
        push!(x,-n)
        push!(y,0)
        for i in 0:1:10
            if (i-5 == 0)
                push!(x, i-5-0.00001)
                push!(x, i-5)
                push!(x, i-5+0.00001)
            else
                push!(x, i-5)
            end
        end
        for i in 0:1:10
            if (i-5 == 0)
                push!(y, 0)
                push!(y, 1)
                push!(y, 0)
            else
                push!(y, 0)
            end
        end
    end
    return [x,y]
end

using CairoMakie
f = Figure()
matrix = kronecker(10)
scatterlines(matrix[1], matrix[2])

##
