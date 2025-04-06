##
function cw_literka_M(width, offset, step)
height = 2
x = Float64[]
y = Float64[]

for i in 1:1:height/step
    push!(x, -width/2+offset)
    push!(y, i)
end
for i in 1:1:width/2/step
    push!(x, -width/2+i+offset)
    push!(y, height/step - height/width/2*i)
end
for i in 1:1:width/2/step
    push!(x, 0+i+offset)
    push!(y, height/width/2*i)
end
for i in 1:1:height/step
    push!(x, width/2+offset)
    push!(y, i)
end

return [x, y]
    
end


matrix =  cw_literka_M(4, 2, 0.01)
##