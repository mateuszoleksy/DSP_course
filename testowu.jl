ramp_wave(t::Real)::Real = 2 * rem(t, 1, RoundNearest)
    sawtooth_wave(t::Real)::Real = -2 * rem(t, 1, RoundNearest)
    triangular_wave(t::Real)::Real = ifelse(mod(t + 1 / 4, 1.0) < 1 / 2, 4mod(t + 1 / 4, 1.0) - 1, -4mod(t + 1 / 4, 1.0) + 3)
    square_wave(t::Real)::Real = ifelse(mod(t, 1) < 0.5, 1, -1)


using CairoMakie
f = lines(range(0,1,1000), ramp_wave.(range(0,1,1000))) 