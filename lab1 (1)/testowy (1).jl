##iodsjfoisdh

t = []
y = []
T = 0.004 
append!(t,0)

for i in 1:10000
    append!(t,T*i)
end

for i in 1:length(t)
    append!(y,sin(200*pi*t[i]))
end

##dgmoidgfo

##quantize(x::Real; L::AbstractVector)::Real = missing

L=[1,2,3,4,5]
x=3.2
y1=[x.>=L]
y2=[x.<=L]