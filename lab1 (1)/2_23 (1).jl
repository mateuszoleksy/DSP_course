##
function heaviside(n)
    x=Float64[]
    y=Float64[]
    if (n < 0)
        push!(x,2*n)
        push!(y,0)
        for i in n*2+1:1:0
            if (i == n)
                push!(x, i-0.000001)
                push!(x, i)
            else
                push!(x, i)
            end
        end
        for i in 2*n:1:0
            if (i >= n)
                push!(y, 1)
            else
                push!(y, 0)
            end
        end
    end
    if (n > 0)
        push!(x,0)
        push!(y,0)
        for i in 1:1:2*n
            if (i == n)
                push!(x, i-0.000001)
                push!(x, i)
            else
                push!(x, i)
            end
        end
        for i in 1:1:2*n+1
            if (i > n)
                push!(y, 1)
            else
                push!(y, 0)
            end
        end
    end
    if (n == 0)
        push!(x,-5)
        push!(y,0)
        for i in -5:1:5
            if (i == 0)
                push!(x, i-0.000001)
                push!(x, i)
                push!(x, i+0.000001)
            else
                push!(x, i)
            end
        end
        for i in -5:1:7
            if (i > n)
                push!(y, 1)
            else
                push!(y, 0)
            end
        end
    end
    return [x,y]
end

using CairoMakie
f = Figure()
matrix = heaviside(0)
scatterlines(matrix[1], matrix[2])
@show matrix[1]
@show matrix[2]
##