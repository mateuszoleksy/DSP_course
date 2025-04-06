module CPS

using LinearAlgebra

author = Dict{Symbol, String}(
    :index => "123456",
    :name  => "Szymon WoĹşniak",
    :email => "szymon.wozniak@agh.edu.pl",
    :group => "0",
)

# SygnaĹy ciÄgĹe
cw_rectangular(t::Real; T=1.0)::Real = missing
cw_triangle(t::Real; T=1.0)::Real = missing
cw_literka_M(t::Real; T=1.0)::Real = missing
cw_literka_U(t::Real; T=1.0)::Real = missing

ramp_wave(t::Real)::Real = missing
sawtooth_wave(t::Real)::Real = missing
triangular_wave(t::Real)::Real = missing
square_wave(t::Real)::Real =  missing
pulse_wave(t::Real, Ď::Real=0.2)::Real = missing
impuse_repeater(g::Function, t1::Real, t2::Real)::Function = missing

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function pulse_wave_bl(t; Ď=0.2, A=1.0, T=1.0, band=20.0)
    missing
end


function impuse_repeater_bl(g::Function, t0::Real, t1::Real, band::Real)::Function
end

function rand_siganl_bl(f1::Real, f2::Real)::Function
    missing
end


# SygnaĹy dyskretne
kronecker(n::Integer)::Real = missing
heaviside(n::Integer)::Real = missing

# Okna
rect(N::Integer)::AbstractVector{<:Real} = missing
triang(N::Integer)::AbstractVector{<:Real} = missing
hanning(N::Integer)::AbstractVector{<:Real} = missing
hamming(N::Integer)::AbstractVector{<:Real} = missing
blackman(N::Integer)::AbstractVector{<:Real} = missing

# Parametry sygnaĹĂłw
mean(x::AbstractVector)::Number = missing
peak2peak(x::AbstractVector)::Real = missing
energy(x::AbstractVector)::Real = missing
power(x::AbstractVector)::Real = missing
rms(x::AbstractVector)::Real = missing

function running_mean(x::AbstractVector, M::Integer)::Vector
    missing
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    missing
end

function running_power(x::AbstractVector, M::Integer)::Vector
    missing
end



# PrĂłbkowanie
function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function=sinc
)
    missing
end

# Kwantyzacja
quantize(L::AbstractVector)::Function = missing
SQNR(N::Integer)::Real = missing
SNR(Psignal, Pnoise)::Real = missing


# Obliczanie DFT
function dtft(f::Real; signal::AbstractVector, fs::Real)
   missing
end

function dft(x::AbstractVector)::Vector
    missing
end

function idft(X::AbstractVector)::Vector
   missing
end

function rdft(x::AbstractVector)::Vector
   missing
end

function irdft(X::AbstractVector, N::Integer)::Vector
   missing
end

function fft_radix2_dit_r(x::AbstractVector)::Vector
   missing
end

function ifft_radix2_dif_r(X::AbstractVector)::Vector
   missing
end

function fft(x::AbstractVector)::Vector
    dft(x) # MoĹźe da rade lepiej?
end

function ifft(X::AbstractVector)::Vector
    idft(X) # MoĹźe da rade lepiej?
end


fftfreq(N::Integer, fs::Real)::Vector = missing
rfftfreq(N::Integer, fs::Real)::Vector = missing
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = missing
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = missing
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = missing

function periodogram(
    x::AbstractVector,
    w::AbstractVector=rect(length(x)),
    L::Integer = 0,
    fs::Real=1.0)::Vector
    missing
end



function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    missing
end


function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    missing
end

function conv(f::Vector, g::Vector)::Vector
    missing
end

function fast_conv(f::Vector, g::Vector)::Vector
    missing
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    missing
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    missing
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    missing
end

function lti_amp(f::Real, b::Vector, a::Vector)::Real
    missing
end

function lti_phase(f::Real, b::Vector, a::Vector)::Real
    missing
end


function firwin_lp_I(order, F0)
    missing
end

function firwin_hp_I(order, F0)
    missing
end

function firwin_bp_I(order, F1, F2)
    missing
end

function firwin_bs_I(order, F1, F2)
    missing
end

function firwin_lp_II(N, F0)
    missing
end

function firwin_bp_II(N, F1, F2)
    missing
end

function firwin_diff(N::Int)
    missing
end

function resample(x::Vector, M::Int, N::Int)
    missing
end



end