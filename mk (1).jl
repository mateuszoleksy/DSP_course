##
# 2.1 ; 2.3 ; 2.7
# 2.10 ; 2.11 ; 2.14 ; 2.15
# 2.22 ; 2.23 ; 2.26
# 3.1 - 3.5
##

## problem 3.1
function mean(x::Vector)
    N = length(x)
    return sum(x)/N
end
x = [1.21 , 23123, 12.32, 2312.321]
mean(x)

## problem 3.2
function peak2peak(x::Vector)
    x1 = maximum(x)
    x2 = minimum(x)
    return abs(x2) + x1
end
x = [1.21 , 23123, 12.32, 2312.321]
peak2peak(x)

## problem 3.3
function energy(x::Vector)
    N = length(x)
    return sum(abs2, x) / N
end
x = [1.21 , 23123, 12.32, 2312.321]
energy(x)

## problem 3.4
function power(x::Vector)
    N = length(x)
    energy_val = sum(abs2, x) / N
    return energy_val^2
end

## problem 3.5
function rms(x::Vector)
    N = length(x)
    rms_val = sum(abs2, x) / N
    return sqrt(rms_val)
end

## problem 2.1
using CairoMakie
A, f, fi, omega , T, n = 2, 25, π/4, 1000, 1/f , 10000
t = collect(0:T/n:T-T/n)
x = A .* sin.(2π.* f .* t .+ fi)

print(x)
plot(t,x)
## problem 2.3
n = 1000
p = 0.25

noise = sqrt(p) * randn(n)
W = noise[1:1000]
## problem 2.7

##