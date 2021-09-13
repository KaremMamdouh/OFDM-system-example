
% clc
% clear all;

%generate the data bits 
N = 5000;
repeat = 0 ; %repeation flag 
col = 2;  %determine number of bits per symbol

% generate random  bits and map them to -1  and 1
%assuming tb and Eo = 1 

data_bits = randi([0 1 ] , N,col);

for k = 1:N
    if(data_bits(k,:) == [0 0])
        data =[data -1-j];
    elseif (data_bits(k,:) == [0 1])
        data(k) = [data -1+j];
    elseif (data_bits(k,:) == [1 0])
        data(k) =[data  1-j];
    elseif (data_bits(k,:) == [1 1])
        data(k) = [data 1+j];
        end
        end

if(repeat > 0) 
    N= N*3;
    data = repelem(data,3,1);
    step = 3; %will be used in descion making
else
    step =1;
end

%channel creation
h = sqrt(0.5) *(normrnd(0,0.5,N,col)+j*normrnd(0,0.5,N,col));

%generation of noise and getting the result  output at
%different No

BER= [];
for S = 1:0.2:20   %signal to noise ratio steps 
    No = 1/S ; 
    %noise generation n 
    n = sqrt(No/(2))*(normrnd(0,0.5,N,col)-j*normrnd(0,0.5,N,col));
    %output formula  
    y = data.*h + n ;
    %there is no need for a correlator since
    %we work with base band signals 
    
    recieved_sig = y./h; %compensate for channel gain
    
    % descion making
for k = 1:col
    for i = 1:step:length(recieved_sig)
 %****** with repeat*******%       
        if(repeat >0)
            if(((recieved_sig(i)>0)&&(recieved_sig(i+1)>0))|| ((recieved_sig(i)>0)&&(recieved_sig(i+2)>0))||((recieved_sig(i+2)>0)&&(recieved_sig(i+1)>0)))
                               
            recieved_sig(i,k) = 1;
            recieved_sig(i+1,k) = 1;
            recieved_sig(i+2,k) = 1;
        else
            recieved_sig(i,k) = -1;
              recieved_sig(i+1,k) = -1;
            recieved_sig(i+2,k) = -1;
        end
        else
%***** no repeat descion****%            
        if(recieved_sig(i,k)>0)
            recieved_sig(i,k) = 1;
        else
            recieved_sig(i,k) = -1;
        end
end
    end
end
%check the BER using built in function symerr
[different,ratio] = symerr(recieved_sig,data,'row-wise');
ratio_symb = mean(ratio);
BER = [BER ratio_symb];
end
SNR = 1:0.2:20;
plot(SNR,BER);
title('BER vs SNR in BPSK');
xlabel('SNR');
ylabel('BER');