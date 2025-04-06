
##
function cw_literka_M(width, height, offset, step)
x = Float64[]
y = Float64[]

for i in 1:1:height/step+1
    push!(x, -width/2/step+offset/step)
    push!(y, i-1)
end
for i in 1:1:width/2/step
    push!(x, -width/2/step+i+offset/step)
    push!(y, height/step - height/(width/2)*i)
end
for i in 1:1:width/2/step
    push!(x, 0+i+offset/step)
    push!(y, height/(width/2)*i)
end
for i in 1:1:height/step+1
    push!(x, width/2/step+offset/step)
    push!(y, i-1)
end

return [x, y]
    
end

using CairoMakie
matrix =  cw_literka_M(4, 4, 0, 0.1)
f = Figure()
f = lines(matrix[1], matrix[2])
f
##