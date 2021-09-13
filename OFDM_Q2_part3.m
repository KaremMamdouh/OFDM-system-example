%%*********QAM*********%%
 clc
 clear all;

%generate the data bits 
N = 5000;
repeat = 0 ; %repeation flag 
col = 1;  %determine number of bits per symbol
SNR = 1:0.2:20;
SNR_db = 10*log10(SNR);
% generate random  bits and map them to -1  and 1
%assuming tb and Eo = 1 

data_bits = randi([0 1 ] , N,col);
data = [];


if(repeat > 0) 
    N= N*3;
    data_bits_  = repelem(data_bits,3,1);

    step = 3; %will be used in descion making
else
    data_bits_ = data_bits;
    step =1;
end


data = qam(data_bits);

%channel creation
h =  (normrnd(0,sqrt(0.5),N/2,1)+j*normrnd(0,sqrt(0.5),N/2,1));

%generation of noise and getting the result  output at
%different No

BERq= [];
recieved_bits = zeros(N/2,2);
for S = 1:0.2:20   %signal to noise ratio steps 
    No = 1/S ; 
    %noise generation n 
    n = normrnd(0,sqrt(No/2),N/2,1)-j*normrnd(0,sqrt(No/2),N/2,1);
    %output formula  
    y = data.*h + n ;
    %there is no need for a correlator since
    %we work with base band signals 
    
    recieved_sig = y./h; %compensate for channel gain
    
    % descion making
recieved_bits = QAM_DEC(recieved_sig);
end
%check the BER using built in function symerr
recieved_bits_ = reshape(recieved_bits',N,1);
[different,ratio] = symerr(recieved_bits_,data_bits_);
BERq = [BERq ratio];
end
SNR = 1:0.2:20;
SNR_db = 10*log10(SNR);

semilogy(SNR_db,BERq);
title('BER vs SNR in 16QAM');
xlabel('SNR dB');
ylabel('BER');