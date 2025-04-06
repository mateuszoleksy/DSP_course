##
function hanning(N)
    return (sin.(range(1,N,N)*pi/N)).^2
end
f = Figure()
lines(range(1,100,100), hanning(100))

##