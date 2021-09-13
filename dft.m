function [ output ] = dft( sig,N )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


if(N > length(sig))
    for i=1:(N-length(sig))
        sig=[sig 0];
    end
end
output=[];
xx=0;
for k=0:(N-1)
    for n=0:N-1
        xx=xx+sig(n+1)* exp(-j*2*pi*n*k/N);
    end
    output=[output xx];
    xx=0;
end

end

