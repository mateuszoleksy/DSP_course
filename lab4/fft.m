x = sin(range(1,1000,1000) * 3.14) 

y = fft(x)

plot(1:1:1000,y,'Color','r')
plot(1:1:1000,x)