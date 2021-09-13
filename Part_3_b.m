%******************** OFDM system simulation ********************** 
clear;clc;close all;
%**************** QAM Flat fading without coding ******************** 

NumberOfBits=128000; % number of bits want to send
NumberOfBitscoded=42000; % number of bits will repeted to send
dataBits=randi([0 1], [1 NumberOfBits]);
dataBitscoded=randi([0 1], [1 NumberOfBitscoded]);
snr_dB=(-4:1:18);
snr_Lin=10.^(snr_dB/10);
BER=0*snr_dB;
BER2=BER;
BERcoded=BER;
BERcoded2=BER;
%**************** coding ******************** 
repeted=repelem(dataBitscoded,3);
coded=[];
for count=0:NumberOfBitscoded*3/63-1
    coded=[coded repeted(count*63+1:(count+1)*63) 0];
end
%**************** interleaver ******************** 
it_data=reshape(dataBits,[8 16 NumberOfBits/128]);
it_datacoded=reshape(coded,[8 16 NumberOfBits/128]);
for k=1:1000
    temp=it_data(:,:,k);
    it_datax(:,:,k)=temp';
    temp=it_datacoded(:,:,k);
    it_datacodedx(:,:,k)=temp';
end
it_data=reshape(it_datax,128,NumberOfBits/128)';
it_datacoded=reshape(it_datacodedx,128,NumberOfBits/128)';
%**************** Mapper ******************** 

% it_data=it_data*2-1;
% it_datacoded=it_datacoded*2-1;
Table=[-3 -1 3 1];
mapped_data=Table(2*it_data(:,1:4:125)+it_data(:,2:4:126)+1)+1i*Table(2*it_data(:,3:4:127)+it_data(:,4:4:128)+1);
mapped_datacoded=Table(2*it_datacoded(:,1:4:125)+it_datacoded(:,2:4:126)+1)+...
    1i*Table(2*it_datacoded(:,3:4:127)+it_datacoded(:,4:4:128)+1);
mapped_data=conj(mapped_data)';
mapped_datacoded=conj(mapped_datacoded)';

%**************** IFFT ******************** 
X_s=ifft(mapped_data,32);
X_scoded=ifft(mapped_datacoded,32);
%**************** cycle prefix ***************** 
X_s=[X_s(end-7:end,:) ; X_s];
X_scoded=[X_scoded(end-7:end,:) ; X_scoded];
%**************** Channel ******************** 
for count=1:length(snr_dB)
    
    noise = sqrt(1/(2*snr_Lin(count))).*(randn(size(X_s))+1i*randn(size(X_s)));
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
     %****** equalizer *********%
    for i = 1:size(rx,2)%executing channel equalization using deconvolution
     y(:,i)=deconv(rx(:,i),flatchanel(:,i));
     y2(:,i)=deconv(rx2(:,i),selective_channel(:,i));
     ycoded(:,i)=deconv(rxcoded(:,i),flatchanel(:,i));
     ycoded2(:,i)=deconv(rxcoded2(:,i),selective_channel(:,i));

    end
    %****** FFT + remove cycle prefix *********%
    x_hat=fft(y(9:end,:),32);  
    x_hat2=fft(y2(9:end,:),32);
    x_hatcoded=fft(ycoded(9:end,:),32);
    x_hatcoded2=fft(ycoded2(9:end,:),32);
    %**************** deMapper ******************** 
    const(:)=[-3-3i, -3-1i, -3+3i, -3+1i,...
                 -1-3i, -1-1i, -1+3i, -1+1i,...
                  3-3i,  3-1i,  3+3i,  3+1i,...
                  1-3i,  1-1i,  1+3i,  1+1i];
    demapped_data=[];
    demapped_data2=[];
    demapped_datacoded=[];
    demapped_datacoded2=[];
for j = 1:32
    for k = 1:1000
    [~, index] = min(abs(x_hat(j,k) - const));
    [~, index2] = min(abs(x_hat2(j,k) - const));
    [~, indexcoded] = min(abs(x_hatcoded(j,k) - const));
    [~, indexcoded2] = min(abs(x_hatcoded2(j,k) - const));
       d(k,:) = dec2bin(index-1,4)-48;
       d2(k,:) = dec2bin(index2-1,4)-48;
       dcoded(k,:) = dec2bin(indexcoded-1,4)-48;
       dcoded2(k,:) = dec2bin(indexcoded2-1,4)-48;
    end
    demapped_data=[demapped_data d];
    demapped_data2=[demapped_data2 d2];
    demapped_datacoded=[demapped_datacoded dcoded];
    demapped_datacoded2=[demapped_datacoded2 dcoded2];
end
    demapped_data=demapped_data';
    demapped_data2=demapped_data2';
    demapped_datacoded=demapped_datacoded';
    demapped_datacoded2=demapped_datacoded2';


%**************** deinterleaver ******************** 
    deit_data=reshape(demapped_data,[16 8 NumberOfBits/128]);%done
    deit_data2=reshape(demapped_data2,[16 8 NumberOfBits/128]);%done
    deit_datacoded=reshape(demapped_datacoded,[16 8 NumberOfBits/128]);%done
    deit_datacoded2=reshape(demapped_datacoded2,[16 8 NumberOfBits/128]);%done

for k=1:1000
    temp=deit_data(:,:,k);
    deit_datax(:,:,k)=temp';
    temp=deit_data2(:,:,k);
    deit_data2x(:,:,k)=temp';
    temp=deit_datacoded(:,:,k);
    deit_datacodedx(:,:,k)=temp';
    temp=deit_datacoded2(:,:,k);
    deit_datacoded2x(:,:,k)=temp';
end
    recieved_bits=deit_datax(:)';
    recieved_bits2=deit_data2x(:)';
    recieved_bitscoded=deit_datacodedx(:)';
    recieved_bitscoded2=deit_datacoded2x(:)';
    
%**************** decoder ******************** 
 recieved_bitsdecoded=[];
 recieved_bitsdecoded2=[];
 for i=0:1999
    recieved_bitsdecoded=[recieved_bitsdecoded round(sum(reshape(recieved_bitscoded(i*64+1:((i+1)*64-1)),3,21))/3)];
    recieved_bitsdecoded2=[recieved_bitsdecoded2 round(sum(reshape(recieved_bitscoded2(i*64+1:((i+1)*64-1)),3,21))/3)];
    
 end

%**************** BER ******************** 
    BER(count)= sum(NumberOfBits-sum(dataBits==recieved_bits))/NumberOfBits;
    BER2(count)= sum(NumberOfBits-sum(dataBits==recieved_bits2))/NumberOfBits;
    BERcoded(count)= sum(NumberOfBitscoded-sum(dataBitscoded==recieved_bitsdecoded))/NumberOfBitscoded;
    BERcoded2(count)= sum(NumberOfBitscoded-sum(dataBitscoded==recieved_bitsdecoded2))/NumberOfBitscoded;
  
end 
 
semilogy(snr_dB,BER,'r*-',snr_dB,BER2,'b*-',snr_dB,(BERcoded),'g*-',snr_dB,(BERcoded2),'c*-','LineWidth',2);

axis([-4 18 10^-2 1]);
title('BER of QPSK  flat fading channel without coding coding');
xlabel('Eb/No (dB)');
ylabel('BER');
grid on;