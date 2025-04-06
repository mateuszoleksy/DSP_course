using LinearAlgebra

mean(x) = sum(x) / length(x)
energy(x) = sum(abs2, x)
power(x) = energy(x) / length(x)
rms(x) = sqrt(power(x))

function poly_from_roots(r)
    p = [1.0 + 0im]
    for i in eachindex(r)
        p = conv(p, [1, -r[i]])
    end
    return p
end

quantize(L) = x -> L[argmin(abs.(-L .+ x))]

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

function lti_phase(b,a,f)
    M = length(b)-1
    K = length(a)-1
    num = sum(cispi(-2*f*m) * b[m+1] for m in 0:M)
    denom = sum(cispi(-2*f*k) * a[k+1] for k in 0:K)
    Hf = num/denom
    return angle(Hf)
end

begin
    function rozwiazanie(;
        zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
        pp::Vector{ComplexF64} = ComplexF64[0.8583027029597319 - 0.4939050783128634im, 0.7335284064927907 + 0.6603284053123151im, 0.8583027029597319 + 0.4939050783128634im, 0.7335284064927907 - 0.6603284053123151im, 0.8164867302004264 - 0.5313204598062244im, 0.7642366912880215 + 0.5986837162256352im, 0.8164867302004264 + 0.5313204598062244im, 0.7642366912880215 - 0.5986837162256352im],
        k::Float64 = 2.2172053767772583e-5,
        F::Vector{Float64} = [0.06, 0.15, 0.2, 0.23, 0.39, 0.41],
    )
        a = poly_from_roots(pp)
        b = poly_from_roots(zz) .* k
        wyniki = [lti_phase(b,a,freq) for freq in F]
        return mean(wyniki)
    end

    out1 = rozwiazanie()
end

begin
    function rozwiazanie(;
        a::Float64 = 0.00026,
        b::Float64 = 0.99,
        x::Vector{Float64} = [0.46899, 0.94554, 0.38483, 0.75908, 0.65304, 0.42565, 0.44475, 0.34957, 0.17759, 0.11466, 0.60465, 0.43936, 0.55578, 0.82076, 0.03235, 0.26533, 0.94833, 0.64894, 0.55192, 0.39964, 0.85519, 0.66929, 0.96504, 0.0008, 0.835, 0.27181, 0.39754, 0.32205, 0.38807, 0.45133, 0.98535, 0.8641, 0.39027, 0.57242, 0.29788, 0.24886, 0.63371, 0.82866, 0.60502, 0.68859, 0.49495, 0.74822, 0.51849, 0.61432, 0.22463, 0.91171, 0.73372, 0.89366, 0.94346, 0.74664, 0.57541, 0.95958, 0.38684, 0.75995, 0.00026, 0.5283, 0.63281, 0.74216, 0.21287, 0.67877, 0.35855, 0.94596, 0.66052, 0.42965, 0.84425],
    )
    N = 8
    L = range(;start=a, stop=b, length=2^N)
        y = quantize(L)
        x_q = y.(x)
        error = x .- x_q
        return rms(error)
    end

    out2 = rozwiazanie()
end


function dft(x)
    N = length(x)
    delta = [cispi(-2 * n /N) for n in 0:N-1]
    [
        sum((
            x[n+1] * delta[(n*k) % N+1] for n in 0:N-1
        )) for k in 0:N-1
    ]
end

function index_to_freq(f, N, fp)
    matrix = f*N/fp
    return Int(round(mod(matrix, N)))+1
end

begin
function rozwiazanie(;
    fp::Int = 1276,
    x::Vector{ComplexF64} = ComplexF64[-0.6 - 0.53im, 0.86 + 0.33im, -0.73 - 0.23im, -0.44 - 0.68im, 0.37 - 1.38im, 0.88 + 0.05im, 1.34 - 0.22im, -0.3 - 0.58im, -0.94 + 1.39im, 0.39 + 1.12im, -0.26 - 0.3im, -0.05 + 0.45im, -0.46 - 1.68im, 0.57 - 0.82im, 0.45 + 0.05im, 0.29 - 1.05im, 0.69 + 1.02im, 1.36 - 0.04im, -1.48 + 0.47im, 0.04 + 0.39im, -0.23 - 1.14im, 1.4 + 1.02im, -0.16 + 0.39im, -0.49 - 0.13im, 0.48 - 0.18im, -0.68 + 0.13im, -1.37 + 0.59im, -0.59 + 0.32im, -0.87 + 0.6im],
    f::Vector{Int} = [-572, 176, 264, 572],
)
    X = dft(x)
    wyniki = [angle(X[index_to_freq(freq, length(x), fp)]) for freq in f]
    return sum(wyniki)
end
    out3 = rozwiazanie()
end


sawtooth_wave(x) = -2 * rem(x, 1.0, RoundNearest)

begin
              
    function rozwiazanie(;
        fp::Float64 = 496.2,
        t1::Float64 = -4.36,
        N::Int = 813,
    )
        x = range(;start=t1, step=(1/fp), length=N)
        y =  0.7*sawtooth_wave.(3.7.*x .- 0.1)
        return mean(y)
    end

    out4 = rozwiazanie()
end

begin
    function rozwiazanie(;
        a::Float64 = -2.6,
        b::Float64 = -1.4,
        x::Vector{Float64} = [-2.48, -2.50968, -2.53936, -2.56905, -2.59873, -1.42841, -1.45809, -1.48777, -1.51746, -1.54714, -1.57682, -1.6065, -1.63618, -1.66587, -1.69555, -1.72523, -1.75491, -1.78459, -1.81428, -1.84396, -1.87364, -1.90332, -1.933, -1.96269, -1.99237, -2.02205, -2.05173, -2.08141, -2.1111, -2.14078, -2.17046, -2.20014, -2.22982, -2.25951, -2.28919, -2.31887, -2.34855, -2.37823, -2.40792, -2.4376, -2.46728, -2.49696, -2.52664, -2.55633, -2.58601, -1.41569, -1.44537, -1.47505, -1.50473, -1.53442, -1.5641, -1.59378, -1.62346, -1.65314, -1.68283, -1.71251, -1.74219, -1.77187, -1.80155, -1.83124, -1.86092, -1.8906, -1.92028, -1.94996, -1.97965, -2.00933, -2.03901, -2.06869, -2.09837, -2.12806, -2.15774, -2.18742, -2.2171, -2.24678, -2.27647, -2.30615, -2.33583, -2.36551, -2.39519, -2.42488, -2.45456, -2.48424, -2.51392, -2.5436, -2.57329, -1.40297, -1.43265, -1.46233, -1.49201, -1.5217, -1.55138, -1.58106, -1.61074, -1.64042, -1.67011, -1.69979],
    )
       N = 9
       L = range(;start=a,stop=b,length=2^N)
       y = quantize(L)
       x_q = y.(x)
       error = x .- x_q
       return rms(error)
    end

    out5 = rozwiazanie()
end

bipolar_square(x) = ifelse(mod(x, 1.0) < 1/2, 1, -1)

    begin
    function rozwiazanie(;
        fp::Float64 = 273.1,
        t1::Float64 = 6.17,
        N::Int = 55,
    )
        x = range(;start=t1, step = (1/fp), length=N)
        y = 3.3 .* bipolar_square.(2.7 .* x .- 2.4)

        return mean(y)
    end

out6 = rozwiazanie()
end

function interpolate(m, s, kernel=sinc)
    return x-> begin
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
        m::Vector{Float64} = [-0.7, -0.6987, -0.6974, -0.6961, -0.6948, -0.6935, -0.6922, -0.6909, -0.6896, -0.6883, -0.687, -0.6857, -0.6844, -0.6831, -0.6818, -0.6805, -0.6792, -0.6779, -0.6766, -0.6753, -0.674, -0.6727, -0.6714, -0.6701, -0.6688, -0.6675, -0.6662, -0.6649, -0.6636, -0.6623, -0.661, -0.6597, -0.6584, -0.6571, -0.6558, -0.6545, -0.6532, -0.6519, -0.6506, -0.6493, -0.648, -0.6467, -0.6454, -0.6441, -0.6428, -0.6415, -0.6402, -0.6389, -0.6376, -0.6363, -0.635, -0.6337, -0.6324, -0.6311, -0.6298, -0.6285, -0.6272, -0.6259, -0.6246, -0.6233, -0.622, -0.6207, -0.6194, -0.6181, -0.6168, -0.6155, -0.6142, -0.6129, -0.6116, -0.6103, -0.609, -0.6077, -0.6064, -0.6051, -0.6038, -0.6025, -0.6012, -0.5999, -0.5986, -0.5973, -0.596, -0.5947, -0.5934, -0.5921, -0.5908, -0.5895, -0.5882, -0.5869, -0.5856, -0.5843],
        s::Vector{Float64} = [0.9165, 0.0479, 0.8957, 0.0641, 0.0749, 0.7044, 0.1083, 0.2622, 0.3433, 0.712, 0.5344, 0.4785, 0.1879, 0.8442, 0.5262, 0.4523, 0.4784, 0.9049, 0.0847, 0.3501, 0.4232, 0.4901, 0.5312, 0.0617, 0.0535, 0.4975, 0.3013, 0.6859, 0.3725, 0.2657, 0.0976, 0.8293, 0.0005, 0.599, 0.937, 0.2805, 0.8521, 0.7385, 0.4645, 0.5344, 0.4411, 0.3273, 0.5718, 0.1282, 0.0584, 0.8397, 0.8474, 0.9712, 0.771, 0.2767, 0.3453, 0.6572, 0.7071, 0.8165, 0.1486, 0.2728, 0.2093, 0.3284, 0.4314, 0.3861, 0.0721, 0.8886, 0.8311, 0.5065, 0.6101, 0.0205, 0.8675, 0.0667, 0.9147, 0.7903, 0.5548, 0.9731, 0.8788, 0.4093, 0.7654, 0.2245, 0.1633, 0.3016, 0.105, 0.8834, 0.7262, 0.882, 0.8665, 0.0007, 0.701, 0.61, 0.1926, 0.8519, 0.69, 0.9246],
        t::Vector{Float64} = [-0.63708, -0.6194, -0.63305, -0.64761, -0.61511],
    )
        y = interpolate(m,s)
        return sum(y.(t))
    end

    out7 = rozwiazanie()
end

begin
function rozwiazanie(;
    a::Float64 = 0.029,
    b::Float64 = 0.98,
    x::Vector{Float64} = [0.60391, 0.95586, 0.20125, 0.04602, 0.07702, 0.73565, 0.83468, 0.07003, 0.7354, 0.81372, 0.82467, 0.67828, 0.72279, 0.83573, 0.42184, 0.19953, 0.35775, 0.27473, 0.6065, 0.84788, 0.31963, 0.37917, 0.643, 0.94691, 0.9553, 0.51686, 0.80469, 0.80766, 0.18759, 0.38219, 0.08759, 0.06497, 0.29981, 0.57455, 0.63507, 0.92241, 0.4596, 0.91797, 0.06605, 0.49719, 0.14322, 0.77556, 0.33664, 0.1476, 0.36395, 0.41187, 0.92945, 0.65771, 0.76112, 0.67985, 0.18495, 0.9815, 0.1183, 0.46507, 0.91572, 0.06659, 0.2457, 0.95056, 0.6229, 0.95903, 0.9342, 0.38223, 0.56734, 0.57739, 0.81512, 0.45378, 0.83719, 0.86254, 0.14285, 0.47133, 0.36332, 0.63342, 0.70516, 0.84559, 0.31345, 0.27463, 0.83197, 0.02881, 0.49899, 0.31653, 0.12998],
)
    N = 3 
    L = range(;start=a, stop=b, length=2^N)
    y = quantize(L)
    x_q = y.(x)
    error = x .- x_q
    return power(error)
end

out8 = rozwiazanie()
end