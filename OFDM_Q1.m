L = 512;
n = 0:2047;
x= sin((2*pi/L)*n);
tic
xf = abs(fft(x));
toc
stem(n,xf);
 tic
 Xf_c = abs(dft(x,length(n)));
 toc
 figure;
 stem(n,Xf_c);