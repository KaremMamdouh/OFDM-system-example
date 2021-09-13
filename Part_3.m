%******************** OFDM system simulation ********************** 
clear;
clc;
close all;
%**************** QPSK Flat fading with  and without coding ******************** 

N_Bits=64000; % number of bits want to send
coded_bits=21000; % repeated bits
dataBits=randi([0 1], [1 N_Bits]); %data bits generated
dataBitscoded=randi([0 1], [1 coded_bits]); %
SNR = (1:1:20);
SNR_dB=10*log10(SNR);
BER_Flat=0*SNR_dB;
%**************** coding ******************** 
bits_repeated=repelem(dataBitscoded,3);
coded=[];
for count=0:coded_bits*3/63-1
    coded=[coded bits_repeated(count*63+1:(count+1)*63) 0];
end
%**************** interleaver ******************** 
inter_data=reshape(dataBits,[8 8 N_Bits/64]);
it_datacoded=reshape(coded,[8 8 N_Bits/64]);
for k=1:1000
    temp=inter_data(:,:,k);
    inter_data(:,:,k)=temp';
    temp=it_datacoded(:,:,k);
    it_datacoded(:,:,k)=temp';
end
inter_data=reshape(inter_data,64,N_Bits/64)';
it_datacoded=reshape(it_datacoded,64,N_Bits/64)';
%**************** Mapper ******************** 
inter_data=inter_data*2-1;
it_datacoded=it_datacoded*2-1;
mapped_data=inter_data(:,1:2:63)+1i*inter_data(:,2:2:64);
mapped_datacoded=it_datacoded(:,1:2:63)+1i*it_datacoded(:,2:2:64);
mapped_data=conj(mapped_data)';
mapped_datacoded=conj(mapped_datacoded)';
%**************** IFFT ******************** 
X_s=ifft(mapped_data,32);
X_scoded=ifft(mapped_datacoded,32);
%**************** cycle prefix ***************** 
X_s=[X_s(end-7:end,:) ; X_s];
X_scoded=[X_scoded(end-7:end,:) ; X_scoded];
%**************** Channel ******************** 
for count=1:length(SNR_dB)
    
    noise = sqrt(1/(2*SNR(count))).*(randn(size(X_s))+1i*randn(size(X_s)));
    flatchanel = sqrt(0.5)*(randn(1,1000)+1i*randn(1,1000));
    selective_channel =sqrt(0.5)*(randn(8,1000)+1i*randn(8,1000));
    X_s_noise=X_s+noise;
    X_s_noisecoded=X_scoded+noise;
    
    for i = 1:size(X_s,2)%convolving the OFDM symbols with the channel and adding AWGN noise
        rx(:,i)=conv(X_s_noise(:,i),flatchanel(:,i));
        rx2(:,i)=conv(X_s_noise(:,i),selective_channel(:,i));
        rxcoded(:,i)=conv(X_s_noisecoded(:,i),flatchanel(:,i));
        rxcoded2(:,i)=conv(X_s_noisecoded(:,i),selective_channel(:,i));

    end
    
%**************** reciever ******************** 
    % companseate for channel effect%
    for i = 1:size(rx,2)%executing channel equalization using deconvolution
     rec_flat(:,i)=deconv(rx(:,i),flatchanel(:,i));
     rec_select(:,i)=deconv(rx2(:,i),selective_channel(:,i));
     rec_flat_coded(:,i)=deconv(rxcoded(:,i),flatchanel(:,i));
     rec_select_coded(:,i)=deconv(rxcoded2(:,i),selective_channel(:,i));

    end
    %****** FFT + remove cycle prefix *********%
    X_F_flat=fft(rec_flat(9:end,:),32);  
    X_F_select=fft(rec_select(9:end,:),32);
    X_F_flat_code=fft(rec_flat_coded(9:end,:),32);
    X_F_select_code=fft(rec_select_coded(9:end,:),32);
    %**************** deMapper ******************** 
 [demapped_data,demapped_select,demapped_flat_code,demapped_select_code] = demapped(X_F_flat,X_F_select,X_F_flat_code,X_F_select_code
%**************** deinterleaver ******************** 
    deint_data=reshape(demapped_data,[8 8 N_Bits/64]);
    deint_data2=reshape(demapped_data2,[8 8 N_Bits/64]);
    deint_datacoded=reshape(demapped_datacoded,[8 8 N_Bits/64]);
    deint_datacoded2=reshape(demapped_datacoded2,[8 8 N_Bits/64]);

for k=1:1000
    temp=deint_data(:,:,k);
    deint_data(:,:,k)=temp';
    temp=deint_data2(:,:,k);
    deint_data2(:,:,k)=temp';
    temp=deint_datacoded(:,:,k);
    deint_datacoded(:,:,k)=temp';
    temp=deint_datacoded2(:,:,k);
    deint_datacoded2(:,:,k)=temp';
end
    recieved_bits=deint_data(:)';
    recieved_bits2=deint_data2(:)';
    recieved_bitscoded=deint_datacoded(:)';
    recieved_bitscoded2=deint_datacoded2(:)';
    
%**************** decoder ******************** 
 recieved_bitsdecoded=[];
 recieved_bitsdecoded2=[];
 for i=0:999
    recieved_bitsdecoded=[recieved_bitsdecoded round(sum(reshape(recieved_bitscoded(i*64+1:((i+1)*64-1)),3,21))/3)];
    recieved_bitsdecoded2=[recieved_bitsdecoded2 round(sum(reshape(recieved_bitscoded2(i*64+1:((i+1)*64-1)),3,21))/3)];
    
 end

%**************** BER ******************** 
    BER_Flat(count)= sum(N_Bits-sum(dataBits==recieved_bits))/N_Bits;
    BER_select(count)= sum(N_Bits-sum(dataBits==recieved_bits2))/N_Bits;
    BER_Flat_code(count)= sum(coded_bits-sum(dataBitscoded==recieved_bitsdecoded))/coded_bits;
    BER_select_code(count)= sum(coded_bits-sum(dataBitscoded==recieved_bitsdecoded2))/coded_bits;
  
end 
 
semilogy(SNR_dB,BER_Flat,SNR_dB,BER_select,SNR_dB,(BER_Flat_code),SNR_dB,(BER_select_code));

title('  BER of QPSK  flat fading channel without coding');
legend('QPSK flat channel   with no repeat','QPSK selective channel   with no repeat','QPSK flat channel   with repeat','QPSK selective channel   with  repeat');
xlabel('SNR dB');
ylabel('BER');
grid on;