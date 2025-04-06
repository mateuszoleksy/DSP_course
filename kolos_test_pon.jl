using LinearAlgebra

mean(x) = sum(x) / length(x)
energy(x) = sum(abs2, x)
power(x) = energy(x)/length(x)
rms(x) = sqrt(power(x))

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

function roots(r)
    H = Matrix(I, length(r)-2, length(r)-2)
    H2 = zeros(length(r)-2)
    H = vcat(H2', H)
    H = hcat(H, -1*reverse(r[2:end]))
    return eigvals(H)
end

begin
    function rozwiazanie(;
        x::Vector{Float64} = [0.17, 1.99, 3.28, 0.51, 3.96, -3.19, -4.27, -2.03, -3.09, -2.3, 4.95, 2.09, -0.2, -0.75, -3.55, -1.72, -0.34, -4.33, 2.15, 0.74, -2.18, 2.07, 2.67, 1.35, -1.76, -0.67, -1.16, 3.34, 3.18, 0.11, -3.49, -2.49, 3.81, 3.78, 2.88, 2.45, 3.61, -2.42, -3.96, -3.02, 4.53, -3.23, 2.37, -4.28, 4.21, -4.12, -3.2, 2.71, -0.12, -1.59, 1.26, -2.12, 2.12, 1.77, 2.57, 2.6, 0.79, 3.07, -0.17, -2.53, -4.42, -2.39, -3.15],
        h::Vector{Float64} = [-4.47, -1.05, -2.91, 0.65, -1.72, 2.0, 2.48, 3.53, -2.18, -2.26, -2.93, -2.01, 2.78, -4.44, -3.62, 4.76],
    )
        wyniki = conv(x,h)
        return mean(wyniki)
    end

    out1 = rozwiazanie()
end

function lti_amp(b,a,f)
    M = length(b)-1
    K = length(a)-1
    num = sum(cispi(-2*m*f) * b[m+1] for m in 0:M)
    denom = sum(cispi(-2*k*f) * a[k+1] for k in 0:K)
    hf = num / denom
    return abs(hf)
end

begin
    function rozwiazanie(;
        b::Vector{Float64} = [0.004445894081944829, 0.017783576327779316, 0.026675364491668976, 0.017783576327779316, 0.004445894081944829],
        a::Vector{Float64} = [1.0, -2.7314228985727755, 3.3881949522875408, -2.1466246391753288, 0.5903324624231749],
        F::Vector{Float64} = [0.06, 0.34, 0.44, 0.47],
    )
        wyniki = [lti_amp(b,a,freq) for freq in F]
        return mean(wyniki)
    end

    out2 = rozwiazanie()
end

begin
    function rozwiazanie(;
        b::Vector{Float64} = [0.011400789339913907, -0.007279142051870437, -0.002889193766565197, 0.0, 0.002889193766565197, 0.00727914205187044, -0.011400789339913906],
        a::Vector{Float64} = [1.0, -1.8026585599670883, 4.308856053853921, -3.6763104388183496, 4.628726419470517, -1.806230871379218, 1.3823464964388523],
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
    out3 = rozwiazanie()
end

kronecker(x) = ifelse(x==0,1,0)
firwin_hp(order, F0) = [Int(kronecker(n)) - 2F0 * sinc(2F0 * n) for n in -order/2:order/2]
hanning(M) = [0.5+0.5cos(2π*n/(2M+1)) for n in -M:M]

begin
    
    function rozwiazanie(;
        order::Int = 94,
        fp::Float64 = 163.0,
        f0::Float64 = 27.71,
        z::Vector{Int} = [68, 23, 23],
    )
        h = firwin_hp(order, f0/fp)
        h = h .* hanning(Int(order/2))
        wyniki = [h[i] for i in z]
        return sum(wyniki)
    end

    out4 = rozwiazanie()
end

function lti_filter(b,a, x)
    M = length(b)
    K = length(a)
    N = length(x)
    sum = zeros(N)
    for n in 1:N
        for m in 0:M-1
            if n-m>0
                sum[n] += b[m+1] * x[n-m]
            end
        end
        for k in 1:K-1
            if n-k>0
                sum[n] -= a[k+1] * sum[n-k]
            end
        end
    end
    return sum
end

begin
    function rozwiazanie(;
        b::Vector{Float64} = [0.03752157066119042, 0.15008628264476168, 0.22512942396714253, 0.15008628264476168, 0.03752157066119042],
        a::Vector{Float64} = [1.0, -1.1813926236104653, 1.3660388577856688, -0.7782356126008878, 0.26718769388569674],
        x::Vector{Float64} = [-0.86, -0.03, -0.14, -0.24, 0.8, 0.55, -0.37, 0.94, 0.52, -0.47, 0.9, -0.08, 0.35],
        L::Int = 26,
    )
        N = length(x)
        x_padded = vcat(x, zeros(Float64, L-N))
        wyniki = lti_filter(b,a,x_padded)
        return rms(wyniki)
    end

    out5 = rozwiazanie()
end

function dft(x)
    N = length(x)
    delta = [cispi(-2*n/N) for n in 0:N-1]
    [
     sum(( x[m+1] * delta[(m*k)%N+1] for m in 0:N-1

     )) for k in 0:N-1
    ]  
end

function freq_to_index(f, N, fp)
    matrix = f*N / fp
    return Int(round(mod(matrix, N)))+1
end

begin
    function rozwiazanie(;
        fp::Int = 1116,
        x::Vector{ComplexF64} = ComplexF64[0.89 - 1.17im, -0.07 - 1.2im, -0.09 + 0.81im, 0.49 + 0.53im, 1.04 + 0.81im, -0.35 - 0.82im, 0.25 - 0.57im, -0.01 + 0.32im, 0.16 - 0.43im, 1.09 - 0.18im, -0.65 + 0.35im, -1.15 - 0.06im, 0.54 + 0.4im, -0.29 - 0.05im, -0.44 + 0.62im, -0.57 - 0.04im, 0.08 + 0.24im, 0.95 + 0.15im, 0.09 + 0.2im, 0.83 + 1.59im, 0.02 - 0.17im, 0.18 - 1.03im, -0.17 - 0.24im, -0.07 - 0.2im, 0.47 + 0.42im, 0.17 + 0.54im, -0.61 + 1.18im, -0.96 + 0.21im, 0.59 - 0.53im, 0.07 - 0.74im, -0.53 - 0.78im],
        f::Vector{Int} = [-540, -468, -288, 0, 252, 396],
    )
    X = dft(x)
    wyniki = [angle(X[freq_to_index(freq, length(x), fp)]) for freq in f]
        return sum(wyniki)
    end

    out6 =rozwiazanie()
end