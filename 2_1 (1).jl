##
function signal(n)
    return 2*sin(2*pi*0.04*n + pi/4)
end

wektor=Float64[]
for i in 1:1:256
push!(wektor, signal(i/1000))
end
|
@show wektor
##  