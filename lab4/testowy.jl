using WAV
using LinearAlgebra

fs = 8e3
t = 0.0:1/fs:prevfloat(1.0)
f = 1e3
y = sin.(2pi * f * t) * 0.1
wavwrite(y, "example.wav", Fs=fs)

y, fs = wavread("example.wav")
y = sin.(2pi * 2f * t) * 0.1
wavappend(y, "example.wav")

y, fs = wavread("example.wav")
wavplay(y, fs)

square_wave(t::Real)::Real = ifelse(mod(t, 1) < 0.5, 1, -1)
using CairoMakie

f = lines(range(1,4000,4000)/1000, square_wave.(range(1,4000,4000)/1000))

function roots(a::AbstractVector)::Vector
    H = Matrix(I, length(a)-2, length(a)-2)
    Z = zeros(length(a) - 2)
    H = vcat(Z', H)
    H = hcat(H, -1 * reverse(a[2:end]))
    @show H
    return eigvals(H)
end
x = [10,15,16,17,18,19,40,50,60]

@show roots(x) 