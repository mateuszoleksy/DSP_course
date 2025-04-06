using LinearAlgebra


quantize(L) = x -> L[argmin(abs.(-L .+ x))]
energy(x) = sum(abs2,x)
power(x) = energy(x) / length(x)
rms(x) = sqrt(power(x))

begin
    function rozwiazanie(;
        a::Float64 = -2.2,
        b::Float64 = -1.1,
        x::Vector{Float64} = [-1.93904, -1.90282, -1.85674, -1.81968, -1.78082, -1.73592, -1.69991, -1.65897, -1.61521, -1.58008, -1.53697, -1.49463, -1.46052, -1.41429, -1.3744, -1.34223, -1.2881, -1.25687, -1.25061, -2.20735, -2.12915, -2.09062, -2.06011, -2.01005, -1.97343, -1.93623, -1.88983, -1.85394, -1.81389, -1.76946, -1.73395, -1.6919, -1.64911, -1.61382, -1.56988, -1.52892, -1.49376, -1.44739, -1.4092, -1.37414, -1.32284, -1.29195, -1.2584, -1.11901, -2.15262, -2.12497, -2.09479, -2.04212, -2.00806, -1.96948, -1.92299, -1.88823, -1.84687, -1.80311, -1.76792, -1.72483, -1.6831, -1.64745, -1.60288, -1.56319],
    )
        N = 2^10
        L = range(;start=a, stop=b, length=N)
        y = quantize(L)
        wyniki = y.(x)
        error = x .- wyniki
        return rms(error)
    end
    out1 = rozwiazanie()
end

begin
    function rozwiazanie(;
        z::Vector{ComplexF64} = ComplexF64[0.5835451923177973 - 0.8120806662658486im, 0.7319786826395969 + 0.6813275336878735im, 0.5835451923177973 + 0.8120806662658486im, 0.7319786826395969 - 0.6813275336878735im, 0.5276490182953228 - 0.8494624850409711im, 0.7673772258842056 + 0.6411959085913295im, 0.5276490182953228 + 0.8494624850409711im, 0.7673772258842056 - 0.6411959085913295im, 1.0 + 0.0im, -1.0 + 0.0im],
        p::Vector{ComplexF64} = ComplexF64[0.6948450536180201 + 0.6999274826683005im, 0.6104754333958993 - 0.7727788955686684im, 0.6948450536180201 - 0.6999274826683005im, 0.6104754333958993 + 0.7727788955686684im, 0.665251017886705 + 0.6843236814847037im, 0.5993315941195554 - 0.738014574026084im, 0.665251017886705 - 0.6843236814847037im, 0.5993315941195554 + 0.738014574026084im, 0.6194566274357605 + 0.6938131908908498im, 0.6194566274357605 - 0.6938131908908498im],
        k::Float64 = 0.003577404602053941,
    )
        p = abs.(p)
        if all(p .< 1)
            return 1.0
        elseif all(p .<= 1) && any(p .== 1)
            return 0.0
        else
            return -1.0
        end
    end

    out2 = rozwiazanie()
end

function dft(x)
    N = length(x)
    delta = [cispi(-2*n/N) for n in 0:N-1]
    [
        sum((   x[m+1] * delta[(m*k) % N+1] for m in 0:N-1

        ))  for k in 0:N-1
    ]
    
end

function freq_to_index(f, N, fp)
    matrix = f*N/fp
    return Int(round(mod(matrix, N)))+1
end

begin
    function rozwiazanie(;
        fp::Int = 624,
        x::Vector{ComplexF64} = ComplexF64[0.63 + 0.26im, -0.3 - 1.14im, 1.79 + 0.83im, -0.89 + 0.52im, -1.04 - 0.72im, -0.51 - 0.29im, -0.55 + 1.19im, -0.84 - 0.13im, 1.12 + 0.66im, -1.1 + 0.18im, 0.24 - 1.76im, 0.16 - 0.64im, -0.28 - 1.05im, -0.46 - 0.45im, -0.67 - 0.03im, 1.03 + 0.11im, 0.14 - 0.4im, -1.26 - 0.28im, 0.89 + 0.31im, 1.05 - 0.3im, 0.45 - 0.36im, 0.32 + 0.44im, -0.34 - 0.37im, 0.01 + 0.73im, -0.58 + 0.46im, -1.57 - 0.29im, -0.4 - 0.31im, 0.23 - 0.35im, 1.27 + 0.57im, -0.64 - 0.51im, -0.43 - 0.31im, -0.23 - 0.9im, 0.58 - 0.98im, 0.33 - 0.64im, -0.31 + 0.2im, -0.03 - 1.09im, -0.19 + 0.99im, -0.42 + 0.13im, 1.39 - 0.81im],
        f::Vector{Int} = [-288, -160, 0, 16, 160],
    )
        X = dft(x)
        phases = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
        return sum(phases)
    end

    out3 = rozwiazanie()
end

firwin_lp(order, F0) = [2F0 * sinc(2F0 * n) for n in -order/2:order/2]
hanning(M) = [0.50 + 0.5*cos(2*pi*n/(2M+1)) for n in -M:M]

begin
    function rozwiazanie(;
        order::Int = 46,
        fp::Float64 = 167.0,
        f0::Float64 = 8.35,
        z::Vector{Int} = [26, 15, 43, 13, 30],
    )
    h = firwin_lp(order, f0/fp)
    h = h .* hanning(Int(order/2))
    wyniki = [h[i] for i in z]
   return sum(wyniki)
    end

    out4 = rozwiazanie()
end


begin
    function rozwiazanie(;
        fp::Int = 448,
        x::Vector{ComplexF64} = ComplexF64[-0.88 - 0.12im, -0.3 - 0.49im, -0.75 - 1.52im, -0.79 - 0.15im, -1.24 - 0.12im, -0.24 + 1.04im, -0.04 - 0.54im, 0.28 - 0.29im, -0.52 - 0.41im, -0.29 + 0.56im, -1.67 - 0.52im, 0.15 - 0.24im, -0.73 + 0.9im, 1.09 + 0.54im, 1.21 + 0.05im, -0.91 - 0.33im, -0.06 - 0.17im, -0.59 + 0.3im, 0.33 - 0.52im, -0.72 + 0.96im, 0.72 - 0.74im, 0.19 - 0.69im, -0.37 - 0.32im, 1.05 - 1.19im, -0.19 - 0.26im, 0.35 - 0.76im, 0.28 - 0.07im, -1.15 + 0.72im],
        f::Vector{Int} = [-160, -144, -80, -32, 80, 144, 224],
    )
    X = dft(x)
    phases = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
    return sum(phases)
    end

    out5 = rozwiazanie()
end

blackmann(M) = [0.42+0.5cos(2π*n/(2M+1))+0.08*cos(4π*n/(2M+1)) for n in -M:M]

begin
    
    function rozwiazanie(;
        fp::Int = 864,
        x::Vector{ComplexF64} = ComplexF64[1.46 + 0.58im, 0.23 + 0.01im, 1.15 - 0.67im, 0.98 - 0.22im, 0.3 + 0.7im, -0.43 + 0.23im, -1.83 - 0.03im, -0.43 - 0.23im, 0.09 - 0.39im, -0.25 - 0.4im, -0.75 + 0.59im, -0.77 + 0.14im, 0.55 - 0.14im, 0.35 - 0.38im, 0.29 - 0.58im, 0.22 - 0.89im, -0.82 - 0.3im, -0.88 + 0.26im, -0.25 + 0.2im, -0.11 + 1.01im, 0.19 - 0.41im, 0.12 + 0.28im, -0.47 - 0.2im, -0.08 + 0.71im, 0.62 + 0.51im, 0.73 + 1.86im, 0.05 - 0.19im, 1.19 - 0.06im, -0.24 + 0.67im, 0.76 + 0.75im, -0.68 - 0.86im, -0.07 - 0.28im, -0.19 - 0.69im, -0.83 - 0.16im, 1.38 + 0.07im, -0.58 - 0.74im],
        f::Vector{Int} = [-384, -336, -240, 264, 432],
    )
    X = dft(x)
    phases = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
    return sum(phases)
    end

    out6 = rozwiazanie()
end

begin
    function rozwiazanie(;
        order::Int = 66,
        fp::Float64 = 192.0,
        f0::Float64 = 23.04,
        z::Vector{Int} = [32, 43, 36, 44, 24, 46],
    )
        h = firwin_lp(order, f0/fp)
        h = h .* blackmann(Int(order/2))
        wyniki = [h[i] for i in z]
        return sum(wyniki)
    end

    out7 = rozwiazanie()
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
        x::Vector{Float64} = [4.71, -2.89, -3.76, -1.46, -2.25, -0.16, 3.02, -1.71, -1.48, -1.56, -4.07, 4.18, -3.97, -0.86, -4.94, 0.26, 1.93, 1.38, 0.44, 3.51, -4.1, 2.63, 2.56, 4.26, -2.53, 0.26, 3.22, 2.9, 3.3, 0.74, 4.37, -1.12, 3.24, 4.77, -4.62, 3.05, 2.19, 1.08, 1.39, -1.45, -2.07, -1.22, 1.79, -3.39, 0.9, -3.81, -4.39, 2.93, -2.52, 3.75, 4.72, -4.89, -2.9, -3.16],
        h::Vector{Float64} = [1.08, 0.37, -2.04, 1.52, -2.88, -4.76, 4.85, -0.42, -1.65, 3.82, 2.56, 0.91, -4.03, -3.59, -3.25, 1.35, 3.94, -1.76, -0.23, 4.34],
    )

    wyniki = conv(x,h)
     return rms(wyniki)
    end

    out8 = rozwiazanie()
end

function interpolate(m,s,kernel::Function=sinc)
    return x-> begin
        sum = 0.0
        delta = m[2] - m[1]
        for i in eachindex(m)
            sum += s[i] * kernel((x - m[i]) / delta )
        end
        return sum
    end
end

begin
    function rozwiazanie(;
        m::Vector{Float64} = [-0.1, -0.0936, -0.0872, -0.0808, -0.0744, -0.068, -0.0616, -0.0552, -0.0488, -0.0424, -0.036, -0.0296, -0.0232, -0.0168, -0.0104, -0.004, 0.0024, 0.0088, 0.0152, 0.0216, 0.028, 0.0344, 0.0408, 0.0472, 0.0536, 0.06, 0.0664, 0.0728, 0.0792, 0.0856, 0.092, 0.0984, 0.1048, 0.1112, 0.1176, 0.124, 0.1304, 0.1368, 0.1432, 0.1496, 0.156, 0.1624, 0.1688, 0.1752, 0.1816, 0.188, 0.1944, 0.2008, 0.2072, 0.2136, 0.22, 0.2264, 0.2328, 0.2392, 0.2456, 0.252, 0.2584, 0.2648, 0.2712, 0.2776, 0.284, 0.2904, 0.2968, 0.3032, 0.3096, 0.316, 0.3224, 0.3288, 0.3352, 0.3416, 0.348, 0.3544, 0.3608, 0.3672, 0.3736, 0.38, 0.3864, 0.3928, 0.3992, 0.4056, 0.412, 0.4184, 0.4248, 0.4312, 0.4376, 0.444, 0.4504, 0.4568, 0.4632, 0.4696, 0.476, 0.4824, 0.4888, 0.4952, 0.5016, 0.508, 0.5144, 0.5208, 0.5272, 0.5336],
        s::Vector{Float64} = [0.8296, 0.2091, 0.971, 0.7679, 0.6197, 0.6184, 0.8686, 0.4627, 0.5325, 0.0221, 0.9989, 0.9498, 0.5218, 0.4512, 0.6459, 0.6648, 0.8564, 0.9205, 0.6908, 0.8649, 0.5851, 0.6414, 0.4494, 0.7905, 0.7084, 0.0374, 0.8872, 0.0436, 0.1427, 0.6865, 0.8375, 0.5243, 0.8889, 0.9137, 0.2664, 0.667, 0.0947, 0.4188, 0.1367, 0.9876, 0.4885, 0.4394, 0.9117, 0.2081, 0.3018, 0.832, 0.6409, 0.8747, 0.4547, 0.1101, 0.2036, 0.5889, 0.7734, 0.9706, 0.4703, 0.2033, 0.1833, 0.1977, 0.2593, 0.924, 0.4564, 0.9378, 0.6854, 0.7679, 0.7975, 0.1627, 0.5981, 0.5746, 0.7448, 0.6904, 0.7923, 0.3591, 0.3733, 0.5592, 0.0289, 0.0667, 0.3027, 0.3272, 0.7326, 0.606, 0.8488, 0.6294, 0.1136, 0.204, 0.0331, 0.8772, 0.2983, 0.172, 0.5312, 0.9459, 0.133, 0.9504, 0.0494, 0.6053, 0.3105, 0.7602, 0.8163, 0.885, 0.6136, 0.7134],
        t::Vector{Float64} = [0.0248, 0.46192, 0.4408, 0.00496, 0.34864, 0.41456, 0.044, 0.53168, 0.4024, 0.17712, 0.2936, 0.50224, 0.3384, 0.13808, 0.41456],
    )
        y = interpolate(m,s)
        return sum(y.(t))
    end

    out9 = rozwiazanie()
end

function lti_filter(b,a,x)
    M = length(b)
    K = length(a)
    N = length(x)
    sum = zeros(N)
    for n in 1:N
        for m in 0:M-1
            if n-m  > 0
                sum[n] += b[m+1] * x[n-m]
            end
        end
        for k in 1:K-1
            if n-k  > 0
                sum[n] -= a[k+1] * sum[n-k]
            end
        end
     end
     return sum
end

begin
    function rozwiazanie(;
        b::Vector{Float64} = [0.00014703394208041414, 0.0, -0.0008822036524824848, 0.0, 0.002205509131206212, 0.0, -0.0029406788416082826, 0.0, 0.002205509131206212, 0.0, -0.0008822036524824848, 0.0, 0.00014703394208041414],
        a::Vector{Float64} = [1.0, -8.337934406429696, 33.0832259783401, -82.42620163610307, 143.42804430133302, -183.46329361908704, 176.7896989620547, -129.2841198089883, 71.22095181816668, -28.841716580428674, 8.159533853997411, -1.4508300482755965, 0.12303590350647985],
        x::Vector{Float64} = [-0.67, 0.3, 0.75, 0.74, -0.41, 0.58, 0.55, -0.27, -0.98, -0.03, 0.52, 0.09, -0.6, 0.21, -0.39, 0.83, -0.89, 0.61, 0.95, 0.8, 0.09, 0.24, 0.52, 0.44, -0.01, 0.95, 0.18, -0.43, -0.22, -0.34, -0.87, -0.23, -0.44, -0.94, -0.16, -0.23, -0.47, -0.29, -0.7, -0.62, 0.98],
        L::Int = 81,
    )
        N = length(x)
        x_padded = vcat(x, zeros(Float64, L-N))
        wyniki = lti_filter(b,a,x_padded)
        return energy(wyniki)
    end

    out10 = rozwiazanie()
end

mean(x) = sum(x)/length(x)

begin
    function rozwiazanie(;
        x::Vector{Float64} = [3.02, 3.48, -2.71, 3.47, -1.17, 4.78, 1.38, 2.2, -4.76, 4.14, -4.22, -0.86, -4.79, 0.11, -5.0, 2.47, 4.62, -2.46, 1.74, -0.87, -1.56, 2.97, 3.7, 3.77, -0.4, 4.12, -4.25, 4.15, -0.46, -1.56, -0.77, -3.0, -4.79, 1.96, 0.24, 0.86, 3.03, 2.85, 0.4, 4.39, 2.89, -2.83, -2.64, 1.34, 1.43, -1.0, 4.22, 0.84, -2.04, 3.78, 4.58, 2.42],
        h::Vector{Float64} = [-1.67, -0.69, -0.18, 3.69, -4.96, -1.75, -4.39, -4.36, 3.56, -3.1, 0.32, -4.76, -2.03, -4.74, 0.82, 3.05, -0.44, 0.34, -2.56, -0.33, 0.54],
    )
        wyniki = conv(x,h)
        return mean(wyniki)
    end

    out11 = rozwiazanie()
end

function roots(r)
    H = Matrix(I, length(r)-2, length(r)-2)
    H = vcat(zeros(length(r)-2)', H)
    H = hcat(H, -reverse(r[2:end]))
    return eigvals(H)
end

begin
    function rozwiazanie(;
        b::Vector{Float64} = [0.034384966351816024, 0.10315489905544807, 0.10315489905544807, 0.034384966351816024],
        a::Vector{Float64} = [1.0, -2.2927784618876466, 2.6459212276108968, -0.9759170014349347],
    )
        p = roots(a)
        p = abs.(p)
        if all(p .< 1)
            return 1.0
        elseif all(p .<= 1) && any(p .== 1)
            return 0.0
        else
            return -1.0
        end
    end

    out12 = rozwiazanie()
end

begin
    function rozwiazanie(;
        x::Vector{Float64} = [-0.83, 0.26, -4.05, -3.47, 4.05, -2.31, -1.1, 1.99, -1.06, -0.79, -4.73, -1.48, 3.12, 1.28, -3.48, -0.69, 2.24, -3.8, -4.3, -0.13, -4.39, 3.57, -2.6, -4.54, 3.45, 0.94, 1.03, -4.94, 0.66, -1.5, 4.12, -0.12, 1.34, 0.84, -1.51, 2.11, -1.53, 2.46, 0.78, -1.1, 4.22, -2.29, 3.27, 0.59, -4.72, 4.86, -3.39, 2.41, 3.52, 0.87, -0.78, 2.26, 4.33, -1.2, -0.94, 2.31, 1.27, -3.91, 0.56, -4.08, -4.75, 2.37, -1.75, -1.35, 2.99, 2.21, 0.65, 0.87, -0.06],
        h::Vector{Float64} = [4.69, -3.6, 3.49, -1.94, 0.24, 4.45, 3.03, 3.57, 4.43, 0.99],
    )
        wyniki = conv(x,h)
        return mean(wyniki)
    end

    out13 = rozwiazanie()
end