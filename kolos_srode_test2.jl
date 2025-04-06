
using LinearAlgebra

firwin_lp(order, F0) = [2F0 * sinc(n*2F0) for n in -order/2:order/2]
hanning(M) = [0.5+0.5*cos(2*pi*n/(2M+1)) for n in -M:M]

begin
    function rozwiazanie(;
        order::Int = 26,
        fp::Float64 = 135.0,
        f0::Float64 = 58.05,
        z::Vector{Int} = [4, 21, 1],
    )
        h = firwin_lp(order, f0/fp)
        h = h .* hanning(Int(order/2))
        hz = [h[i] for i in z]
        return sum(hz)
    end

    out1 = rozwiazanie()
end

function roots(a)
    H = Matrix(I, length(a)-2, length(a)-2)
    H = vcat(zeros(Float64, length(a)-2)', H)
    H = hcat(H, -reverse(a[2:end]))
    return eigvals(H)
end

begin    
    function rozwiazanie(;
        b::Vector{Float64} = [5.939790904914056e-5, 0.00023759163619656225, 0.0003563874542948434, 0.00023759163619656225, 5.939790904914056e-5],
        a::Vector{Float64} = [1.0, -6.639026700511004, 13.608552706804907, -11.332260038315336, 3.36431562304722],
    )
        p = roots(a)
        p = abs.(p)
        if all(p .< 1.0)
            return 1.0
        elseif all(p .<= 1.0) && any(p .== 1.0)
            return 0.0
        else
            return -1.0
        end

    end

    out2 = rozwiazanie()
end

quantize(L) = x -> L[argmin(abs.(-L .+ x))]
energy(x) = sum(abs2, x)
power(x) = energy(x) / length(x)
rms(x) = sqrt(power(x))

begin
    function rozwiazanie(;
        a::Float64 = -0.89,
        b::Float64 = 1.5,
        x::Vector{Float64} = [0.54, 0.6384, 0.73679, 0.83519, 0.93358, 1.03198, 1.13037, 1.22877, 1.32717, 1.42556, -0.87604, -0.77765, -0.67925, -0.58086, -0.48246, -0.38406, -0.28567, -0.18727, -0.08888, 0.00952, 0.10791, 0.20631, 0.30471, 0.4031, 0.5015, 0.59989, 0.69829, 0.79668, 0.89508, 0.99348, 1.09187, 1.19027, 1.28866, 1.38706, 1.48545, -0.81615, -0.71775, -0.61936, -0.52096, -0.42257, -0.32417, -0.22578, -0.12738, -0.02898, 0.06941, 0.16781, 0.2662, 0.3646, 0.46299, 0.56139, 0.65979, 0.75818, 0.85658, 0.95497, 1.05337, 1.15176, 1.25016, 1.34856, 1.44695, -0.85465, -0.75626, -0.65786, -0.55947, -0.46107, -0.36267, -0.26428, -0.16588, -0.06749, 0.03091, 0.1293, 0.2277, 0.3261, 0.42449, 0.52289, 0.62128, 0.71968, 0.81807, 0.91647, 1.01487, 1.11326, 1.21166, 1.31005, 1.40845, -0.89316, -0.79476, -0.69636],
    )
        N = 2^8
        L = range(;start=a, stop=b, length=N)
        y = quantize(L)
        x_quantized = y.(x)
        error = x .- x_quantized
        return rms(error)
    end

    out3 = rozwiazanie()
end

function lti_amp(b,a,f)
    M = length(b)
    K = length(a)
    num = sum([b[m+1]*cispi(-2*m*f) for m in 0:M-1])
    denom = sum([a[k+1]*cispi(-2*k*f) for k in 0:K-1])
    hf = num/denom
    return abs(hf)
end

mean(x) = sum(x) / length(x)

begin
    
    function rozwiazanie(;
        b::Vector{Float64} = [0.1910246184379355, -0.9551230921896775, 1.910246184379355, -1.910246184379355, 0.9551230921896775, -0.1910246184379355],
        a::Vector{Float64} = [1.0, -1.9190105608014512, 2.0469666180566546, -0.9665161989736166, 0.25152894005410054, 0.07123452787189105],
        F::Vector{Float64} = [0.04, 0.11, 0.14, 0.45, 0.46, 0.46],
    )
        wyniki = [lti_amp(b,a,freq) for freq in F]
        return mean(wyniki)
    end
     out4 = rozwiazanie()
end