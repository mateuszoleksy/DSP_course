x = 1:1:1000
y = sin(x)


y1 = fft(y,length(y))

plot(x,y1)