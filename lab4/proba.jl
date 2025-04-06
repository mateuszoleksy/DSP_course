
##
# Funkcja do obliczania FFT
function fftm(x)
    N = length(x)
    if N <= 1
        return x
    else
        even = fft(x[1:2:end])
        odd = fft(x[2:2:end])
        factor = exp.(-2im * π * (0:N÷2-1) / N)
        return vcat(even + factor .* odd, even - factor .* odd)
    end
end

# Testowanie implementacji
using CairoMakie
frequencies = fftm(sin.(range(1,100,100)))
plot(range(1,100,100), frequencies)
##