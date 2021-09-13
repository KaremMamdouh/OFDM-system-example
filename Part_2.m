%******************** BPSK,QPSK,16-PSK,8-QAM Initialization ********************** 
clear;clc;close all;

Framelength=51200; % number of symbols you want to send

snr_Lin=1:1:20;
snr_dB=10*log10(snr_Lin);

Es_N0_dB_16  = snr_dB + 10*log10(log2(16));%***16-QAM
Es_N0_dB_4 = snr_dB + 10*log10(log2(4));%****QPSK

BER=zeros(1,length(snr_dB));
BER4=zeros(1,length(snr_dB));
BER16q=zeros(1,length(snr_dB));

%******************** for Reptation ********************* 
BERrp=zeros(1,length(snr_dB));
BER4rp=zeros(1,length(snr_dB));
BER16qrp=zeros(1,length(snr_dB));

%******************** START BER CALCULATION ********************* 
loops=1;  % Number of simulation loops to increase the accuracy of ber

for i=1:length(snr_dB)
    BERtemp=zeros(1,loops);
    BERtemp4=zeros(1,loops);
    BERtemp16q=zeros(1,loops);
    
    BERtempre=zeros(1,loops);
    BERtemp4re=zeros(1,loops);
    BERtemp16qre=zeros(1,loops);
    for ii=1:loops 
     
%******************** Data generation ********************************   
 
    data=randi([0 1],1,Framelength);
    dataq=randi([0 3],1,Framelength);
    data16=randi([0 15],1,Framelength);
    
%****************** Repetied Data generation *************************
    datare=repelem(data,3);
    bin=dec2bin(dataq,2)';
    dataqre=bin2dec((reshape(repelem(bin(:),3),2,Framelength*3)'))';
    bin=dec2bin(data16,4)';
    data16re=bin2dec(reshape(repelem(bin(:),3),4,Framelength*3)')';   
    
%******************** Mapper ***********************   
    
    %********Bpsk***************%
    theta=2*pi*data/2;%modulation turning data to angle
    tx=-1*real(exp(1i*theta)); %converting phasor into real and image
    %**Repetation**%
    txre=2*datare-1;
    
    %********QPSK**********%
    datagreyqbits=fliplr(de2bi(dataq))';
    datagreyqbits=2*(datagreyqbits)-1;
    % two bpsk signal one in-phase and the other quadrature phase
    si=datagreyqbits(1,:);%***in-phase part we send half of data    
    sq=datagreyqbits(2,:);%***Quadrature part we send next half of data
    txq=si+1j*sq;
    
    %**Repetation**%
    datagreyqbits=fliplr(de2bi(dataqre))';
    datagreyqbits=2*(datagreyqbits)-1;
    % two bpsk signal one in-phase and the other quadrature phase
    sire=datagreyqbits(1,:);%***in-phase part we send half of data    
    sqre=datagreyqbits(2,:);%***Quadrature part we send next half of data
    txqre=sire+1j*sqre;
    
    %********16-QAM******%
    %***Constellation Table******% 
    table2(:) = [-3-3i, -3-1i, -3+3i, -3+1i,...
                 -1-3i, -1-1i, -1+3i, -1+1i,...
                  3-3i,  3-1i,  3+3i,  3+1i,...
                  1-3i,  1-1i,  1+3i,  1+1i];
    %***mapping data to constellation table values    
    tx16q = table2(data16+1);
    tx16qre = table2(data16re+1);
%************ AWGN Channel ***************
	%********Creating Noise signal with defined power******%
    nc=random('norm',0,sqrt(1/(snr_Lin(i)*2)),[1 Framelength]);
    ns=random('norm',0,sqrt(1/(snr_Lin(i)*2)),[1 Framelength]);
    noisesignal=nc + 1i*ns;
    noisesignalq=sqrt(0.5*(10^((10*log10(log2(4))-Es_N0_dB_4(i))/10)))*...
                            (randn(1,Framelength)+1i*randn(1,Framelength));  
    tx16q_power = 10*log10(sum(abs(tx16q(:)).^2)/length(tx16q(:)));
    noisesignal16q=sqrt((0.5*(10.^((tx16q_power-Es_N0_dB_16(i))/10))))*(randn(1,Framelength)+1i*randn(1,Framelength));

       %**Repetation**%
    noisesignalre=sqrt(1/(snr_Lin(i)*2))*(randn(1,Framelength*3)+1i*randn(1,Framelength*3));
    noisesignalqre=sqrt(0.5*(10^((10*log10(log2(4))-Es_N0_dB_4(i))/10)))*...
                            (randn(1,Framelength*3)+1i*randn(1,Framelength*3));
    tx16q_power = 10*log10(sum(abs(tx16qre(:)).^2)/length(tx16qre(:)));
    noisesignal16qre=sqrt((0.5*(10.^((tx16q_power-Es_N0_dB_16(i))/10))))*(randn(1,Framelength*3)+1i*randn(1,Framelength*3));
    
	%********Creating Rayleigh flat fading channel effect ******%
    hr=random('norm',0,1/sqrt(2),[1 Framelength]);
    hi=random('norm',0,1/sqrt(2),[1 Framelength]);
    %h=1;
    h=hr + 1i*hi;
        %**Repetation**%
    hr=random('norm',0,1/sqrt(2),[1 Framelength*3]);
    hi=random('norm',0,1/sqrt(2),[1 Framelength*3]);
    hre=hr + 1i*hi;
    
    %****adding noise and Channel effect to transmitted signal*****%
    rx=tx.*h+noisesignal;
    rxq=txq.*h+noisesignalq;
    rx16q=tx16q.*h+noisesignal16q;

        %**Repetation**%
    rxre=txre.*hre+noisesignalre;
    rxqre=txqre.*hre+noisesignalqre;
    rx16qre=tx16qre.*hre+noisesignal16qre;


%******************** Demapper **********************
    %******BPSK*******************%
     rx=rx./h;
     rx=rx(:);
     datademodbp=(((rx) > 0)+0)';
     
             %**Repetation**%
     rxre=rxre./hre;
     rxre=rxre(:);
     datademodbpre=(((rxre) > 0)+0)';
          
    %**********QPSK******************%
    rxq=rxq./h;
    si_=sign(real(rxq));                                
    sq_=sign(imag(rxq));
    ber1=(Framelength-sum(si==si_))/Framelength; %In-phase BER calculation
    ber2=(Framelength-sum(sq==sq_))/Framelength; %Quadrature BER calculation    
    
            %**Repetation**%
    rxqre=rxqre./hre;
    datademodqre=(sign(real(rxqre))+1)+(1+sign(imag(rxqre)))/2;
    
     %*******16-QAM*******************%
     %****Constellation Table****%
     const2=table2(:);
     rx16q=rx16q./h;   
     rx16q=rx16q(:);
for N = 1:length(rx16q)
     %compute the minimum distance from each symbol in constellation
     %minimum distance defines to us the location of symbol in
     %constellation
     [~, index] = min(abs(rx16q(N) - const2));
       datademod16q(N) = index-1;
end
   
        %**Repetation**%
     rx16qre=rx16qre./hre;   
     rx16qre=rx16qre(:);
for N = 1:length(rx16qre)
          [~, index] = min(abs(rx16qre(N) - const2));
       datademod16qre(N) = index-1;
end

%******************** Derepeter ******************
    
    %******BPSK*******************%
    datademodbpre=round(sum(reshape(datademodbpre,3,Framelength))/3);
    %******QPSK*******************%
    bin=dec2bin(datademodqre,2)';
    bin=round(sum(reshape(bin(:)-48,3,Framelength*2))/3);
    datademodqre=bin2dec(char(reshape(bin+48,2,Framelength)'))';
    %******16-QAM*******************%
    bin=dec2bin(datademod16qre,4)';
    bin=round(sum(reshape(bin(:)-48,3,Framelength*4))/3);
    datademod16qre=bin2dec(char(reshape(bin+48,4,Framelength)'))'; 
    
%******************** Bit Error Rate (BER) ******************
    BERtemp(ii)= sum(Framelength-sum(data==datademodbp))/length(data);
    BERtemp4(ii)=mean([ber1 ber2]);
    BERtemp16q(ii)=  sum(Framelength-sum(de2bi(data16)==de2bi(datademod16q)))/(length(data16)*4);
    
    BERtempre(ii)   = sum(Framelength-sum(data==datademodbpre))/length(data);
    BERtemp4re(ii)  = sum(Framelength-sum(de2bi(dataq)==de2bi(datademodqre)))/(length(data16)*2);
    BERtemp16qre(ii)= sum(Framelength-sum(de2bi(data16)==de2bi(datademod16qre)))/(length(data16)*4);

    end     
  BER(i)=sum(BERtemp)/loops; %average ber for loops
  BER4(i)=sum(BERtemp4)/loops;
  BER16q(i)=sum(BERtemp16q)/loops;
  
  BERre(i)=sum(BERtempre)/loops; %average ber for loops
  BER4re(i)=sum(BERtemp4re)/loops;
  BER16qre(i)=sum(BERtemp16qre)/loops;

end

%****************Theoretical BER Calculations*************************%
TheoryBER=0.5*erfc(sqrt(10.^(snr_dB/10)));
TheoryBER4=0.5*erfc(sqrt(10.^(snr_dB/10)));
TheoryBER16q=3/8*erfc(sqrt(snr_Lin/2.5));

FadingBER=berfading(snr_dB,'psk',2,1);
FadingBER4=berfading(snr_dB,'psk',4,1);
FadingBER16q=berfading(snr_dB,'QAM',16,1);
% **********************  BER Plot ***************************% 
figure(1);
% semilogy(snr_dB,BER,'r*-',snr_dB,BER4,'c+-',snr_dB,BER16q,'b+-','LineWidth',2);
semilogy(snr_dB,BER4);

% hold on;
% semilogy(snr_dB,TheoryBER,'k.-',snr_dB,TheoryBER4,'y--',snr_dB,TheoryBER16q,'m--','LineWidth',2);
% axis([-4 16 10^-4 1]);
% title('AWGN BER Theoretical/Simulation for BPSK,QPSK,16-QAM without coding');
% legend('BPSK','QPSK','16-QAM','BPSK Theoretical','QPSK Theoretical','16-QAM Theoretical');
% xlabel('Eb/No (dB)');
% ylabel('BER');
% grid on;
% hold off
%%
figure(2);
semilogy(snr_dB,BER,'r*-',snr_dB,BER4,'c+-',snr_dB,BER16q,'b+-','LineWidth',2);
hold on;
semilogy(snr_dB,FadingBER,'k.-',snr_dB,FadingBER4,'y--',snr_dB,FadingBER16q,'m--','LineWidth',2);
axis([-4 16 10^-3 1]);
title('Raylrigh BER Theoretical/Simulation for BPSK,QPSK,16-QAM without coding');
legend('BPSK','QPSK','16-QAM','BPSK Theoretical Fading','QPSK Theoretical Fading','16-QAM Theoretical Fading');
xlabel('Eb/No (dB)');
ylabel('BER');
grid on;
hold off


figure(3);
semilogy(snr_dB,BER16qre);
title('BER vs SNR in 16QAM');
legend('16QAM with  repeat');

 xlabel('SNR dB');
 ylabel('BER');
 hold on
%  semilogy(snr_dB,BERre,'k*-',snr_dB,BER4re,'g+-',snr_dB,BER16qre,'m+-','LineWidth',2);
% axis([-4 16 10^-3 1]);
% title('Raylrigh BER Simulation for BPSK,QPSK,16-QAM coding Vs without coding');
% legend('BPSK','QPSK','16-QAM','BPSK Theoretical','QPSK Theoretical','16-QAM Theoretical');
% legend('BPSK no coding','QPSK no coding','16-QAM no coding','BPSK with coding','QPSK with coding','16-QAM with coding');
% xlabel('Eb/No (dB)');
% ylabel('BER');
% grid on;
hold off
