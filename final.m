%% Q1

clear;
N=2048;
xi = randn(1,N);
tic;
Xi_dft = dft(xi,N);
t_dft = toc

tic;
Xi_fft = fft(xi,N);
t_fft=toc
%% Q2
% a) BPSK No repetition code

no_rand_bits = 5000;
bk = randi(2,1,no_rand_bits)-1;
%% 
% 


Eb = 1;
rep = 1;

xk = BPSKGen(bk,rep,Eb);
no_rand_bits = length(xk);
%generating the random channel

h = rayleighGen(no_rand_bits);

%generating the random Noise
Eb_No = -50:10:50;
s=1;
BER = [];
for k = Eb_No
    n = noiseGen(k, no_rand_bits);
    %Compute the received symbol vector
    yk = xk.*h+n;
    %Compensate for the channel gain at the receiver
    bk_hat = BPSK_Decision(yk,h,no_rand_bits);
    bit_error = ber_calc(xk,bk_hat,rep);
    BER(s) = bit_error;
    s=s+1;
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
title('BPSK over AWGN channel');
%% 
% b) BPSK with repetition code
%%

Eb = 1;
rep = 3;

xk = BPSKGen(bk,rep,Eb);
no_rand_bits = length(xk);
%generating the random channel

h = rayleighGen(no_rand_bits);

%generating the random Noise
Eb_No = -50:10:50;
s=1;
BER = [];
for k = Eb_No
    n = noiseGen(k, no_rand_bits);
    %Compute the received symbol vector
    yk = xk.*h+n;
    %Compensate for the channel gain at the receiver
    bk_hat = BPSK_Decision(yk,h,no_rand_bits);
    bit_error = ber_calc(xk,bk_hat,rep);
    BER(s) = bit_error;
    s=s+1;
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
title('BPSK over AWGN channel with repetition');
%% 
% c) QPSK no repetition
%%

Eb = 1;
rep = 1;

xk = QPSKGen(bk,rep);
no_rand_bits = length(xk)*rep;
%generating the random channel

h = rayleighGen(no_rand_bits);

%generating the random Noise
Eb_No = -50:10:50;
s=1;
BER = [];
for k = Eb_No
    n = noiseGen(k, no_rand_bits);
    %Compute the received symbol vector
    yk = xk.*h+n;
    %Compensate for the channel gain at the receiver
    bk_hat = QPSK_Decision(yk,h);
    bit_error = ber_calc(xk,bk_hat,rep);
    BER(s) = bit_error;
    s=s+1;
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
title('QPSK over AWGN channel with no repetition');
%% 
% d) QPSK with repition 
%%

Eb = 1;
rep = 3;

xk = QPSKGen(bk,rep);
no_rand_bits = length(xk);
%generating the random channel
h = rayleighGen(no_rand_bits);

%generating the random Noise
Eb_No = -50:10:50;
s=1;
BER = [];
for k = Eb_No
    n = noiseGen(k, no_rand_bits);
    %Compute the received symbol vector
    yk = xk.*h+n;
    %Compensate for the channel gain at the receiver
    bk_hat = QPSK_Decision(yk,h);
    bit_error = ber_calc(xk,bk_hat,rep);
    BER(s) = bit_error;
    s=s+1;
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
title('QPSK over AWGN channel with repetition');
%% 
% e) 16-QAM with no repition 


Eb = 1;
rep = 1;

xk = QAMGen(bk,rep);
no_rand_bits = length(xk);
%generating the random channel
h = rayleighGen(no_rand_bits);

%generating the random Noise
Eb_No = -50:10:50;
s=1;
BER = [];
for k = Eb_No
    n = noiseGen(k, no_rand_bits);
    %Compute the received symbol vector
    yk = xk.*h +n;
    %Compensate for the channel gain at the receiver
    bk_hat = QAM_Decision(yk,h);
    bit_error = ber_calc(xk,bk_hat,rep);
    BER(s) = bit_error;
    s=s+1;
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
title('QAM over AWGN channel with no repetition');
%% 
% f) 16-QAM with repition 


Eb = 1;
rep = 3;

xk = QAMGen(bk,rep);
no_rand_bits = length(xk);
%generating the random channel
h = rayleighGen(no_rand_bits);

%generating the random Noise
Eb_No = -50:10:50;
s=1;
BER = [];
for k = Eb_No
    n = noiseGen(k, no_rand_bits);
    %Compute the received symbol vector
    yk = xk.*h +n;
    %Compensate for the channel gain at the receiver
    bk_hat = QAM_Decision(yk,h);
    bit_error = ber_calc(xk,bk_hat,rep);
    BER(s) = bit_error;
    s=s+1;
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
title('QAM over AWGN channel with repetition');
%% Q3
% a) QPSK Simulation
%%
technique_len=64; no_of_bits = 64; no_of_symbols = 100; rep = 1;
channel_type =1 ;
xk = bitsGeneration(technique_len,no_of_bits,no_of_symbols,rep);
xk_inter = interleaver(xk,8,8, no_of_symbols);
xk_mapped = mapper(xk_inter,no_of_symbols)
xk_ifft = ifft(xk_mapped);
xk_cp = CyclicExt(xk_ifft);
Eb_No=1:50
for EbNo = Eb_No
    [yk,h] = channel(xk_cp,EbNo,30);
%     plot(abs(fftshift(yk(:,:,1))))
%     hold on;
%     plot(abs(fftshift(h(:,:,1))))
%     
%     hold off ;
    mod_scheme = 1;
    [yk_hat] = reciever(yk,h,mod_scheme)
    yk_hat = interleaver(yk_hat,8,8,no_of_symbols);
    bit_error = ber_calc(xk,yk_hat,rep)
    BER(EbNo)=bit_error 
end
semilogy(Eb_No,BER);grid;xlabel('SNR');ylabel('BER');
%% Functions
%%
%%%%%%%% Q1 Functions %%%%%%%%%%%
%%
function [X] = dft(x,N)
%dft Summary of this function goes here
%   Detailed explanation goes here
k=0:N-1;
w= (1/N) * 2* pi * k;
d = 1;
X =0;
for n = 0:length(x)-1
    e = x(d)*exp (-1i.*n.*w);
    X = X + e;
    d=d+1;
end
end

%%%%%%%%%%%%%% Q2 Functions %%%%%%%%%%%
function x_BPSK = BPSKGen(bk,rep,Eb)
bk_nrz = bk*2 - 1;
x_BPSK = bk_nrz*sqrt(Eb);
x_BPSK = repelem(x_BPSK,rep);
end

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

function x_QAM = QAMGen(bk,rep)
QAM16_Symbols = [ -3-3i -3-1i -3+3i -3+1i -1-3i -1-1i -1+3i -1+1i 1-3i 1-1i 1+3i 1+1i 3-3i 3-1i 3+3i 3+1i];
bk = repelem(bk,rep);
x_QAM = zeros(1,length(bk)/4);
k = 1;
for g = 1:length(bk)/4
    idx = (bk(k)*8+bk(k+1)*4+bk(k+2)*2+bk(k+3))+1;
    x_QAM(g) = QAM16_Symbols(idx);
    k =k+4;
end
end

function [h] = rayleighGen(len)
%generating the random channel
hr = randn(1,len)*sqrt(0.5);
hi = 1j*randn(1,len)*sqrt(0.5);
h = hr + hi;
end
% 
% function [h] = rayleighGen(len)
% %generating the random channel
% taps =1 ;
% h = [(randn(1,taps)+1i*randn(1,taps))*sqrt(0.5) zeros(1,len-taps)].' ; % generate channel taps and pad with zeros
% end
% 

function [n] = noiseGen(Eb_No, len)
Eb_No_lin = 10.^(0.1*Eb_No);
nc = randn(1,len)*sqrt(0.5/Eb_No_lin);
ns = 1j*randn(1,len)*sqrt(0.5/Eb_No_lin);
n=nc+ns;
end

function bk_hat = BPSK_Decision(yk,h,no_rand_bits)
yk = yk./h;
bk_hat = ones(1,no_rand_bits);
for v = 1:no_rand_bits
    if real(yk(v)) <= 0
        bk_hat(v) = -1;
    end
end
end

function bk_hat = QPSK_Decision(yk,h)
yk = yk ./ h;
no_rand_bits = length(yk);
bk_hat = ones(1,no_rand_bits);
for v = 1:no_rand_bits
    if real(yk(v)) <= 0 && imag(yk(v)) <= 0
        bk_hat(v) = -1 - 1i;
    elseif  real(yk(v)) <= 0 && imag(yk(v)) >= 0
        bk_hat(v) = -1 + 1i;
    elseif  real(yk(v)) >= 0 && imag(yk(v)) <= 0
        bk_hat(v) = 1 - 1i;
    elseif  real(yk(v)) >= 0 && imag(yk(v)) >= 0
        bk_hat(v) = 1 + 1i;
    end
end
end

function bk_hat = QAM_Decision(yk,h)
yk = yk ./ h;
no_rand_bits = length(yk);
bk_hat = ones(1,no_rand_bits);
for v = 1:no_rand_bits
    if real(yk(v)) <= 0 && real(yk(v)) >= -2 &&  imag(yk(v)) <= 0 && imag(yk(v)) >= -2
        bk_hat(v) = -1 - 1i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) <= -2 &&  imag(yk(v)) <= 0 && imag(yk(v)) >= -2
        bk_hat(v) = -3 - 1i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) >= 2 &&  imag(yk(v)) <= 0 && imag(yk(v)) >= -2
        bk_hat(v) = 3 - 1i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) <= 2 &&  imag(yk(v)) <= 0 && imag(yk(v)) >= -2
        bk_hat(v) = 1 - 1i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) >= -2 &&  imag(yk(v)) >= 0 && imag(yk(v)) <= 2
        bk_hat(v) = -1 + 1i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) <= -2 &&  imag(yk(v)) >= 0 && imag(yk(v)) <= 2
        bk_hat(v) = -3 + 1i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) >= 2 &&  imag(yk(v)) >= 0 && imag(yk(v)) <= 2
        bk_hat(v) = 3 + 1i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) <= 2 &&  imag(yk(v)) >= 0 && imag(yk(v)) <= 2
        bk_hat(v) = -1 + 1i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) >= -2 &&  imag(yk(v)) <= 0 && imag(yk(v)) <= -2
        bk_hat(v) = -1 - 3i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) <= -2 &&  imag(yk(v)) <= 0 && imag(yk(v)) <= -2
        bk_hat(v) = -3 - 3i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) >= 2 &&  imag(yk(v)) <= 0 && imag(yk(v)) <= -2
        bk_hat(v) = 3 - 3i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) <= 2 &&  imag(yk(v)) <= 0 && imag(yk(v)) <= -2
        bk_hat(v) = 1 - 3i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) >= -2 &&  imag(yk(v)) >= 0 && imag(yk(v)) >= 2
        bk_hat(v) = -1 + 3i;
    elseif  real(yk(v)) <= 0 && real(yk(v)) <= -2 &&  imag(yk(v)) >= 0 && imag(yk(v)) >= 2
        bk_hat(v) = -3 + 3i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) >= 2 &&  imag(yk(v)) >= 0 && imag(yk(v)) >= 2
        bk_hat(v) = 3 + 3i;
    elseif  real(yk(v)) >= 0 && real(yk(v)) <= 2 &&  imag(yk(v)) >= 0 && imag(yk(v)) >= 2
        bk_hat(v) = 1 + 3i;
    end
end
end

function bit_error = ber_calc(xk,bk_hat,rep)
bit_error = 0;
count = 0;
no_rand_bits = length(xk);
for v = 1:rep:no_rand_bits
    for s = 0 : rep-1
        if xk(v+s) ~= bk_hat(v+s)
            count = count + 1;
        end
    end
    if count >= rep / 2
        bit_error = bit_error + 1;
    end
    count = 0;
end
bit_error = bit_error / no_rand_bits;
end

%%%%%%%%%%%%%%%%%%%%%% Q3 Functions %%%%%%%%%%%%%%%%%%%%%%
function xk = bitsGeneration(technique_len,no_of_bits,no_of_symbols,rep)
rng('default');
rng(1);
xk = zeros(1,technique_len,no_of_symbols);
for r = 1:no_of_symbols
    bk = randi(2,1,no_of_bits)-1;
    bk = repelem(bk,rep);
    if length(bk)<64
        bk  = [bk  zeros(1,technique_len-length(bk))];
    else
        bk = bk(1:technique_len);
    end
    xk(:,:,r) = bk;
end
end

function xk_inter=interleaver(xk,sz1,sz2,no_of_symbols)
for r = 1: no_of_symbols
    x=xk(:,:,r);
    xs(:,:,r) = reshape(x,sz1,sz2);
    xk_inter(:,:,r) = reshape(xs(:,:,r).' ,1,[]);
end
end

function xk_mapped = mapper(xk_inter,no_of_symbols)
mod_technique = size(xk_inter);
if mod_technique(2) == 64 %QPSK
    QPSK_Symbols = [ -1-1i -1+1i 1-1i 1+1i ];
    xk_mapped = zeros(1,32,no_of_symbols);
    for r= 1:no_of_symbols
        k = 1;
        for g = 1:32
            idx = (xk_inter(k)*2+xk_inter(k+1))+1;
            xk_mapped(:,g,r) = QPSK_Symbols(idx);
            k =k+2;
        end
    end
elseif mod_technique(2) == 128 %16-QAM
    QAM16_Symbols = [ -3-3i -3-1i -3+3i -3+1i -1-3i -1-1i -1+3i -1+1i 1-3i 1-1i 1+3i 1+1i 3-3i 3-1i 3+3i 3+1i];
    xk_mapped = zeros(1,32,no_of_symbols);
    for r= 1:no_of_symbols
        k = 1;
        for g = 1:32
            idx = (xk_inter(k)*8+xk_inter(k+1)*4+xk_inter(k+2)*2+xk_inter(k+3))+1;
            xk_mapped(:,g,r) = QAM16_Symbols(idx);
            k =k+4;
        end
    end
end
end

function [xk_cp] = CyclicExt(xk_ifft)
%CYCLICEXT Summary of this function goes here
%   Detailed explanation goes here
sz= size(xk_ifft);
size_of_sig= sz(2);
xk_cp = zeros(sz(1),sz(2)*1.25,sz(3));
for r = 1:sz(3)
    xk_cp(:,:,r) = [xk_ifft(:,(3*size_of_sig/4+1):size_of_sig,r) xk_ifft(:,:,r)];
end
end

function [yk,h] = channel(xk_cp,EbNo,taps)
%CHANNEL Summary of this function goes here
%   Detailed explanation goes here
rng('default');
rng(1);
sz= size(xk_cp);
for r = 1: sz(3)
    len= sz(2);
    h(:,:,r) = [randn(1,taps)+j*randn(1,taps) zeros(1,sz(2)-taps+1)].' ; % generate channel taps and pad with zeros
    Eb_No_lin = 10.^(0.1*EbNo);
    nc = randn(1,len)*sqrt(0.5/Eb_No_lin);
    ns = 1j*randn(1,len)*sqrt(0.5/Eb_No_lin);
    n = nc + ns ;
    xn = xk_cp(:,:,r) + n;
    yk(:,:,r) = conv(h(:,:,r),xn).';
end
end



function [yk_hat] = reciever(yk,h,mod_scheme)
%RECIEVER Summary of this function goes here
%   Detailed explanation goes here
sz= size(yk);
yk_cp = [];
for r = 1:sz(3)
    yk_rec(:,:,r) = deconv(yk(:,:,r), h(:,:,r));% deconv to remove the channel effect
    sz= size(yk_rec);
    size_of_sig= sz(2);
    yk_cp(:,:,r) = yk_rec(:,size_of_sig/5+1:size_of_sig,r);%remove the cyclic prefix
end

size(yk_cp);
yk_fft = fft(yk_cp); % revert to freq domain
sz = size(yk_fft);

for r = 1:sz(3)
    if mod_scheme ==1
        yk_hat = zeros(1,sz(2)*2,sz(3));
        l = 1;
        for v = 1:sz(2)%32
            if real(yk_fft(1,v,r)) <= 0 && imag(yk_fft(1,v,r)) <= 0
                yk_hat(1,l,r) = 0;
                yk_hat(1,l+1,r) = 0;
            elseif  real(yk_fft(1,v,r)) <= 0 && imag(yk_fft(1,v,r)) >= 0
                yk_hat(1,l,r) = 0;
                yk_hat(1,l+1,r) = 1;
            elseif  real(yk_fft(1,v,r)) >= 0 && imag(yk_fft(1,v,r)) <= 0
                yk_hat(1,l,r) = 1;
                yk_hat(1,l+1,r) = 0;
            elseif  real(yk_fft(1,v,r)) >= 0 && imag(yk_fft(1,v,r)) >= 0
                yk_hat(l) = 1;
                yk_hat(l+1) = 1;
            end
            l = l + 2;
        end
    elseif mod_scheme ==2
        yk_hat = zeros(1,sz(2)*2,sz(3));
        
        
    end
end

end