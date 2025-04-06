using LinearAlgebra


function dft(x)
    N = length(x)
    delta = [cispi(-2*n/N) for n in 0:N-1]
    [
        sum(( x[m+1] * delta[(m*k)%N+1] for m in 0:N-1

        )) for k in 0:N-1
    ]
    
end

function freq_to_index(f, N, fp)
    matrix = f*N/fp
return Int(round(mod(matrix, N)))+1
    
end

begin
    function rozwiazanie(;
        fp::Int = 572,
        x::Vector{ComplexF64} = ComplexF64[1.73 - 0.24im, -0.25 - 0.31im, 0.02 - 0.32im, 0.65 + 0.92im, -1.14 + 0.05im, 0.3 - 0.05im, -0.42 + 0.75im, 0.27 - 1.38im, 1.0 - 0.18im, 0.87 - 0.38im, -1.07 - 0.52im, 0.17 + 0.24im, -0.48 + 0.39im, -0.72 - 1.22im, -0.78 + 0.67im, 1.06 - 0.01im, -0.13 - 0.51im, 1.16 - 0.24im, -0.72 - 0.77im, -0.73 + 1.18im, -0.75 + 0.07im, -0.59 - 0.96im, -1.0 + 0.05im, -0.24 + 0.65im, -0.0 + 0.6im, -0.58 - 0.81im, -0.52 + 0.72im, -0.98 + 0.66im, -0.17 + 1.1im, 1.03 - 0.05im, -0.15 + 0.37im, 1.14 - 0.52im, -0.07 + 0.57im, -1.25 - 0.28im, -0.15 + 1.99im, 0.08 + 0.85im, -0.14 + 0.05im, -0.62 + 0.17im, -0.59 - 0.37im, 0.78 + 0.44im, 0.4 - 0.1im, -0.67 - 0.31im, 0.47 - 0.11im, -1.31 - 0.26im],
        f::Vector{Int} = [-247, -195, -117, -39, 26, 117],
    )
        X = dft(x)
        wyniki = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
        return sum(wyniki)
    end

    out1 = rozwiazanie()
end

square_wave(x) = ifelse(mod(x, 1.0) < 1/2, 1, -1)
energy(x) = sum(abs2, x)
power(x) = energy(x) / length(x)
rms(x) = sqrt(power(x))

 begin
    function rozwiazanie(;
        fp::Float64 = 232.19,
        t1::Float64 = -1.28,
        N::Int = 481,
    )
        g = square_wave
        x = range(;start=t1, step=1/fp, length=N)
        y = 1.1*g.(0.3*x.-4.4)
        return energy(y)
    end

    out2 = rozwiazanie()
end

quantize(L) = x-> L[argmin(abs.( -L .+ x))]

begin
    function rozwiazanie(;
        a::Float64 = -5.0,
        b::Float64 = 1.1,
        x::Vector{Float64} = [-2.896, -2.49831, -2.10062, -1.70294, -1.30525, -0.90756, -0.50987, -0.11218, 0.2855, 0.68319, 1.08088, -4.92143, -4.52375, -4.12606, -3.72837, -3.33068, -2.93299, -2.53531, -2.13762, -1.73993, -1.34224, -0.94455, -0.54687, -0.14918, 0.24851, 0.6462, 1.04388, -4.95843, -4.56074, -4.16305, -3.76536, -3.36768, -2.96999, -2.5723, -2.17461, -1.77692, -1.37924, -0.98155, -0.58386, -0.18617, 0.21151, 0.6092, 1.00689, -4.99542, -4.59773, -4.20005, -3.80236, -3.40467, -3.00698, -2.60929, -2.21161, -1.81392, -1.41623, -1.01854, -0.62086, -0.22317, 0.17452, 0.57221, 0.9699],
    )
        L = range(;start=a, stop=b, length=2^8)
        y = quantize(L)
        x_quantized = y.(x)
        error = x .- x_quantized
        return rms(error)
    end

    out3 = rozwiazanie()
end

function interpolate(m,s,kernel::Function=sinc)
    return x -> begin
        sum = 0.0
        delta = m[2] - m[1]
        for i in eachindex(m)
            sum += s[i] * kernel((x - m[i]) / delta)
        end
        return sum
    end
end

begin
    function rozwiazanie(;
        m::Vector{Float64} = [-4.6, -4.5993, -4.5986, -4.5979, -4.5972, -4.5965, -4.5958, -4.5951, -4.5944, -4.5937, -4.593, -4.5923, -4.5916, -4.5909, -4.5902, -4.5895, -4.5888, -4.5881, -4.5874, -4.5867, -4.586, -4.5853, -4.5846, -4.5839, -4.5832, -4.5825, -4.5818, -4.5811, -4.5804, -4.5797, -4.579, -4.5783, -4.5776, -4.5769, -4.5762, -4.5755, -4.5748, -4.5741, -4.5734, -4.5727, -4.572, -4.5713, -4.5706, -4.5699, -4.5692, -4.5685, -4.5678, -4.5671, -4.5664, -4.5657, -4.565, -4.5643, -4.5636, -4.5629],
        s::Vector{Float64} = [0.5958, 0.3847, 0.6169, 0.7813, 0.6945, 0.0898, 0.4739, 0.8492, 0.931, 0.2699, 0.3124, 0.8272, 0.6107, 0.1566, 0.6121, 0.14, 0.1262, 0.9958, 0.9448, 0.6351, 0.7042, 0.9129, 0.745, 0.0746, 0.3708, 0.6316, 0.5439, 0.0429, 0.4309, 0.0559, 0.2413, 0.7904, 0.3413, 0.919, 0.5045, 0.6059, 0.6965, 0.5052, 0.6718, 0.7854, 0.3825, 0.199, 0.3226, 0.5655, 0.9585, 0.2263, 0.2492, 0.8288, 0.7006, 0.5951, 0.3114, 0.52, 0.6104, 0.6906],
        t::Vector{Float64} = [-4.59615, -4.56304, -4.59993, -4.5706, -4.58264, -4.6, -4.5867, -4.5692, -4.5762, -4.57088],
    )
        y = interpolate(m,s)
        wyniki = y.(t)
        return sum(wyniki)
    end

    out4 = rozwiazanie()
end

begin 
    function rozwiazanie(;
        m::Vector{Float64} = [0.8, 0.8054, 0.8108, 0.8162, 0.8216, 0.827, 0.8324, 0.8378, 0.8432, 0.8486, 0.854, 0.8594, 0.8648, 0.8702, 0.8756, 0.881, 0.8864, 0.8918, 0.8972, 0.9026, 0.908, 0.9134, 0.9188, 0.9242, 0.9296, 0.935, 0.9404, 0.9458, 0.9512, 0.9566, 0.962, 0.9674, 0.9728, 0.9782, 0.9836, 0.989, 0.9944, 0.9998, 1.0052, 1.0106, 1.016, 1.0214, 1.0268, 1.0322, 1.0376, 1.043, 1.0484, 1.0538, 1.0592, 1.0646, 1.07, 1.0754, 1.0808, 1.0862, 1.0916, 1.097, 1.1024, 1.1078, 1.1132, 1.1186, 1.124, 1.1294, 1.1348, 1.1402, 1.1456],
        s::Vector{Float64} = [0.5879, 0.0624, 0.338, 0.3487, 0.966, 0.9881, 0.4618, 0.9414, 0.0557, 0.2742, 0.8242, 0.7884, 0.2368, 0.2029, 0.3251, 0.7451, 0.444, 0.2463, 0.3765, 0.9548, 0.9549, 0.1737, 0.4429, 0.1661, 0.4207, 0.4321, 0.1567, 0.2096, 0.1712, 0.0889, 0.5548, 0.0676, 0.0647, 0.9079, 0.2596, 0.9944, 0.7446, 0.5548, 0.6366, 0.3511, 0.526, 0.0669, 0.4058, 0.2677, 0.289, 0.1582, 0.6312, 0.8197, 0.6836, 0.8976, 0.6642, 0.5623, 0.6249, 0.8976, 0.0862, 0.501, 0.2032, 0.4405, 0.1057, 0.8131, 0.4032, 0.6804, 0.2005, 0.4042, 0.0417],
        t::Vector{Float64} = [0.95498, 0.83078, 0.84428, 0.91448, 1.0727, 0.86912, 0.8702, 1.11968, 0.95876, 0.95336, 0.82106, 0.86156, 1.02248, 1.01978, 1.05596, 1.03058],
    )
        y = interpolate(m,s)
        wyniki = y.(t)
        return sum(wyniki)
    end

    out5 = rozwiazanie()
end

begin
    function rozwiazanie(;
        fp::Int = 1974,
        x::Vector{ComplexF64} = ComplexF64[0.62 + 0.16im, -0.63 - 0.06im, 0.32 - 0.54im, 1.14 - 0.22im, -0.62 - 0.35im, -0.66 + 0.74im, 0.81 - 0.75im, 1.48 - 0.77im, -0.14 - 0.31im, -0.4 + 0.02im, -0.44 - 0.77im, -0.57 - 1.32im, 0.1 - 0.85im, -0.04 - 0.05im, -0.82 - 0.56im, -0.59 + 0.82im, 0.6 - 2.45im, 0.32 - 0.66im, -0.75 - 0.36im, 0.08 - 0.21im, -0.2 - 0.04im, 0.87 - 1.17im, -1.17 + 0.3im, -0.35 + 0.55im, -0.48 + 0.32im, 1.2 + 0.19im, 0.11 - 0.06im, 0.35 + 1.59im, -0.62 + 0.97im, 1.01 + 0.38im, -1.11 - 0.31im, -0.78 + 0.14im, -0.84 + 0.87im, 0.01 - 0.06im, 0.01 - 0.92im, -1.04 - 0.31im, 0.4 + 0.8im, 1.07 - 1.6im, 0.24 + 0.67im, 0.08 - 0.79im, -0.17 - 0.35im, -0.63 + 0.67im, 0.22 + 1.13im, -0.0 + 0.37im, -1.45 - 0.16im, -0.4 - 1.15im, -1.31 + 0.83im],
        f::Vector{Int} = [-924, -42, 336, 756, 798, 924],
    )
        X = dft(x)
        wyniki = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
        return sum(wyniki)
    end

    out6 = rozwiazanie()
end
function conv(f,g)
    M = length(f)
    K = length(g)
    sum = zeros(eltype(f), M+K-1)
    for m in 1:M
        for k in 1:K
            sum[m+k-1] += f[m] * g[k]
        end
    end
    return sum
end

begin
    function rozwiazanie(;
        x::Vector{Float64} = [-0.12, 4.53, 3.22, 3.48, 4.0, -0.54, 3.47, 3.66, -1.38, -2.73, 4.78, 1.55, -3.62, 2.92, -4.84, -2.54, 1.89, 4.22, -2.02, -2.07, 4.81, -0.29, -4.16, 3.0, -2.22, 1.35, -3.81, 2.07, -3.29, 1.94, -3.69, -4.05, 0.76, -0.21, -2.85, 4.3, -0.67, -2.34, 3.47, 0.88, 0.15, -3.78, 2.77, -1.28, 1.56, -4.56, -4.64, -1.53, -0.43, 0.4, 1.96, 1.14, 2.75, 4.59, 4.97, -2.98, -4.87, -0.97, 0.96, 1.25, 4.91, -3.17, -2.06, 0.83, -3.04, -3.22, -1.41, 1.12],
        h::Vector{Float64} = [-1.06, 2.51, -1.3, -4.18, 0.18, 0.16, 0.2, -3.01, -1.41, 1.55, -0.9, 1.94, 2.67, 4.81, 3.63, -4.06, -1.12, 1.55, 3.31, 4.61, -3.67],
    )
        y = conv(x,h)
        return power(y)
    end

    out7 = rozwiazanie()
end


begin
    function rozwiazanie(;
        fp::Int = 375,
        x::Vector{ComplexF64} = ComplexF64[-0.3 + 0.31im, 0.17 + 0.07im, -0.04 - 0.23im, 0.17 - 0.27im, 0.53 + 1.71im, 0.07 - 1.14im, 0.1 - 0.89im, 0.14 + 1.12im, -0.71 + 0.03im, -2.1 + 1.28im, -0.41 - 1.38im, -0.59 + 0.61im, -0.55 + 0.61im, -1.35 + 0.92im, -0.17 - 0.49im, -0.2 - 0.61im, -0.34 - 0.6im, -0.2 + 0.09im, -0.74 + 0.49im, -0.29 - 0.12im, -0.09 - 0.59im, -0.9 - 0.73im, -0.98 - 0.79im, 0.14 + 0.17im, 1.15 + 1.51im],
        f::Vector{Int} = [-150, -75, 45, 60, 75, 165],
    )
        X = dft(x)
        wyniki = [abs(X[freq_to_index(freq, length(x), fp)]) / length(x) for freq in f]
        return sum(wyniki)
    end

    out8= rozwiazanie()
end

begin
    function rozwiazanie(;
        fp::Int = 1536,
        x::Vector{ComplexF64} = ComplexF64[0.22 + 0.17im, 0.03 - 0.72im, 1.03 + 0.95im, 0.68 + 0.78im, 0.02 + 0.37im, -0.21 - 0.57im, -0.29 + 0.6im, -1.07 + 1.03im, 0.33 - 0.34im, 0.38 - 0.69im, -0.56 + 0.18im, 0.49 + 0.02im, -0.06 - 1.23im, -1.2 + 0.07im, -0.13 + 0.31im, 0.2 - 1.01im, 1.26 + 1.22im, 0.5 + 0.44im, 0.34 + 0.83im, 1.21 - 0.16im, -1.41 + 0.43im, 0.47 + 1.26im, 0.42 - 0.57im, 0.18 + 0.58im, -1.11 + 0.92im, 0.07 + 0.36im, 0.85 - 0.61im, 0.04 + 0.51im, 1.25 + 0.66im, -1.17 + 0.32im, -0.53 - 0.61im, -0.34 + 0.62im, -0.49 - 0.56im, 0.06 + 0.12im, 0.03 + 0.64im, -0.66 + 0.44im, -0.25 + 0.04im, 0.91 - 0.22im, 0.41 - 0.15im, -0.74 - 0.69im, 0.8 - 0.5im, 0.95 - 0.51im, -0.01 - 0.29im, -0.27 + 0.53im, 0.52 - 1.06im, -0.34 + 0.72im, -0.04 + 0.08im, -0.02 - 0.51im],
        f::Vector{Int} = [-416, 64, 416, 640, 768],
    )
    X = dft(x)
    wyniki = [abs(X[freq_to_index(freq, length(x), fp)]) / length(x) for freq in f]
    return sum(wyniki)
    end

    out9 = rozwiazanie()
end

begin
    function rozwiazanie(;
        fp::Int = 1564,
        x::Vector{ComplexF64} = ComplexF64[2.09 - 0.33im, 0.02 - 0.41im, -0.79 + 0.54im, -1.22 + 0.44im, 0.97 + 0.43im, 0.35 - 1.04im, 0.19 - 0.45im, 0.85 - 0.3im, 1.23 - 0.29im, 0.73 - 0.68im, -0.12 + 0.61im, 0.43 + 1.0im, -0.58 - 0.48im, 0.26 - 0.87im, -0.79 - 1.19im, -0.26 - 0.18im, 0.35 + 0.15im, 0.68 - 0.33im, -0.32 - 0.45im, -0.12 - 0.43im, 0.57 - 0.79im, 0.19 - 1.05im, 0.69 - 1.46im, -0.09 - 0.08im, -0.15 + 0.44im, -0.26 - 0.85im, -0.01 + 0.88im, 1.23 - 1.02im, 0.45 - 1.16im, -0.15 - 0.28im, -1.05 + 1.35im, -0.43 + 0.27im, -1.06 - 1.14im, 0.52 - 0.95im],
        f::Vector{Int} = [-690, -552, -460, -230, 46, 690],
    )
    X = dft(x)
    wyniki = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
    return sum(wyniki)
    end
    out10 = rozwiazanie()
end

begin
    function rozwiazanie(;
        fp::Int = 1540,
        x::Vector{ComplexF64} = ComplexF64[0.45 + 0.88im, 0.56 + 0.42im, -0.89 + 1.56im, 0.69 + 0.83im, -0.6 + 2.1im, 0.08 + 0.12im, -0.43 + 1.25im, 0.08 - 0.23im, 1.3 - 0.95im, 0.4 - 0.08im, -0.17 - 0.55im, 1.12 - 0.27im, -0.69 + 0.44im, -0.1 + 0.47im, -0.17 - 0.7im, 1.07 - 0.66im, -0.32 - 0.07im, 0.45 + 0.1im, 0.22 + 0.32im, 0.06 + 0.19im, -1.42 - 0.13im, -1.0 - 0.12im, -0.99 - 0.45im, -0.15 + 0.29im, 0.19 + 0.22im, 0.21 + 0.12im, -0.86 + 1.0im, 1.05 - 0.86im, 0.47 + 0.3im, -0.47 - 1.03im, -1.4 - 0.02im, -0.37 + 0.39im, -1.5 + 0.43im, 0.3 + 0.25im, 0.63 + 0.04im, 1.59 + 0.63im, 0.03 - 0.51im, -0.48 + 0.12im, 1.2 - 0.77im, 0.81 + 0.13im, 0.76 - 0.65im, 0.36 - 0.11im, 1.14 - 0.8im, -0.92 + 0.04im],
        f::Vector{Int} = [-700, -455, -315, 175, 315, 770],
    )
    X = dft(x)
    wyniki = [abs(X[freq_to_index(freq, length(x), fp)]) / length(x) for freq in f]
    return sum(wyniki)
    end

    out11 = rozwiazanie()
end

function lti_filter(b,a,x)
    M = length(b)
    K = length(a)
    N = length(x)
    sum = zeros(N)
    for n in 1:N
        for m in 0:M-1
            if n-m > 0
                sum[n] += b[m+1] * x[n-m]
            end
        end
        for k in 1:K-1
            if n-k > 0
                sum[n] -= a[k+1] * sum[n-k]
            end
        end
    end
    return sum
end

mean(x) = sum(x)/length(x)

begin
    function rozwiazanie(;
        b::Vector{Float64} = [0.34405884576897705, -0.9956140623053321, 0.9956140623053321, -0.344058845768977],
        a::Vector{Float64} = [1.0, -0.9647025660886663, 0.6005193136649185, -0.1141239363950345],
        x::Vector{Float64} = [-1.0, -0.72, 0.28, -0.28, 0.06, 0.51, 0.21, -0.35, 0.03, 0.19, 0.61, -0.57, 0.46, 0.77, -0.73, 0.62, 0.71, -0.41, 0.89, -0.15, 0.45, 0.09, -0.02, -0.25, -0.97, 0.98, -0.06, 0.31, -0.6, 0.59, -0.67, 0.33, 0.79],
        L::Int = 54,
    )
        N = length(x)
        x_padded = vcat(x, zeros(Float64, L-N))
        y = lti_filter(b,a,x_padded)
        return mean(y)
    end

    out12 = rozwiazanie()
end