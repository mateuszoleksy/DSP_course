##
    function cw_literka_M(t::Int64, T::Int64)
        step = 0.01
        return [[((ones(1,floor(Int, T/step)).-1).+floor(Int, t/step)),
        (range(0,floor(Int, T/2/step+1),floor(Int, T/2/step+1))).+floor(Int, t/step),
        (range(floor(Int, T/2/step+1),floor(Int, T/step),floor(Int, T/2/step+1))).+floor(Int, t/step), 
        ((ones(1,floor(Int, T/step)).-1).+floor(Int, t/step).+floor(Int, T/step))], [ range(1,floor(Int, T/step),floor(Int, T/step)), 
        (((range(0,floor(Int, T/2/step), floor(Int, T/2/step+1)).+floor(Int, t/step))).*(-1)).+floor(Int, t/step).+floor(Int, T/step), 
        (((range(floor(Int, T/2/step),floor(Int, T/step), floor(Int, T/2/step+1)).+floor(Int, t/step))*(1)).-floor(Int, t/step)), 
        range(floor(Int, T/step),1,floor(Int, T/step))]
        ]
    end

    using CairoMakie
    f, ax, l1 = lines(vec(cw_literka_M(2,4)[1][1]), vec(cw_literka_M(2,4)[2][1]))
    l2 = lines!(ax, vec(cw_literka_M(2,4)[1][2]), vec(cw_literka_M(2,4)[2][2]))
    l3 = lines!(ax, vec(cw_literka_M(2,4)[1][3]), vec(cw_literka_M(2,4)[2][3]))
    l4 = lines!(ax, vec(cw_literka_M(2,4)[1][4]), vec(cw_literka_M(2,4)[2][4]))
    f
##