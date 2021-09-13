function x_QPSK = QPSKGen(bk,rep)
bk = repelem(bk,rep);
QPSK_Symbols = [-1-i -1+i 1-i 1+i];
x_QPSK = zeros(1,length(bk)/2);
k=1;
for b = 1:length(bk)/2
    idx = (bk(k)*2+bk(k+1))+1;
    x_QPSK(b) = QPSK_Symbols(idx) ;
end
end
