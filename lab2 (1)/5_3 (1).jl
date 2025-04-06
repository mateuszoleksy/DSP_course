function SNR(Psignal, Pnoise)
    return 20*log10(Psignal/Pnoise)
end

@show SNR(5, 0.25)

