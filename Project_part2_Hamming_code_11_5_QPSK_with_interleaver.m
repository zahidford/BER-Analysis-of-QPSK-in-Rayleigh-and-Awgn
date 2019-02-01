clc
close all;
clear all;

%Simulation of coded QPSK system in base band (AWGN + Rayleigh fading channel)
%Channel coding: single error correcting code Hamming (15,11)
%With interleaver

%Coherence and signal time calculation
n=15; %length of the codeword
k=11; %length of the message
t=1; %Number of corrected bits 
fc=1e10; %operating frequency (Hz)
c=3e8; %light speed
lambda=c/fc; %wavelength of the carrier
v=6e4/3600; %Speed (m/s) == 60 km/h    
fm=v/lambda; %Maximum Doppler frequency
Tcoh=(9/16/pi)*(1/fm); %Coherence time
Rb=1e6; %Bit rate (bit/s)
Tb=1/Rb; %Bit duration (s/bit)
Ts=2*Tb; %Signal time (QPSK)
d=round(Tcoh/(Tb*(k/n)))+1; %Depth of the interleaver (greater than Tcoh/Tb-coded)

N=1936000; %Total number of bits of the data stream 
EBN0dB=0:1:20; %Normalized SNR per bit in dB
EBN0=10.^(EBN0dB/10); %Normalized SNR per bit 
ECN0dB = EBN0dB + 10*log10(k/n); %Normalized SNR per coded bit in dB
ECN0=10.^(ECN0dB/10); %Normalized SNR per coded bit
N0=10.^(-ECN0dB/10); %Noise spectral density

I=eye(k);
P1=[1 0 1 1 1 0 0 0 1 1 1]';
P2=[1 1 0 1 1 0 1 1 0 0 1]';
P3=[1 1 1 0 1 1 0 1 1 0 0]';
P4=[1 1 1 1 0 1 1 0 0 1 0]';
P=[P1 P2 P3 P4];
G=[P I];
H=[eye(n-k) transpose(P)];

%Theoretical coded QPSK BER with Rayleign fading and interleaving
Pc=(1/2)*(1-sqrt(ECN0./(1+ECN0)));
aux=0; %Auxiliar variable
theory_coded_interleaved_BER=0;

for p=(t+1):1:n
aux=(1/n)*(p*(factorial(n)/(factorial(p)*factorial(n-p)))*(Pc.^p).*(1-Pc).^(n-p));
theory_coded_interleaved_BER=theory_coded_interleaved_BER+aux;
end

theory_coded_BER=(1/2)*(1-sqrt(ECN0./(1+ECN0))); %Theoretical coded QPSK BER with Rayleign fading

for x=1:length(EBN0dB) 
    data=round(rand(1,N)); %Random data stream

    i=data(1:2:end); %In-phase bits
    q=data(2:2:end); %Quadrature bits
    
    ui=reshape(i,k,N/k/2)';        %in-phase 11-bits message
    uq=reshape(q,k,N/k/2)';        %Quadrature 11-bits message
    
    ci=mod(ui*G,2);              %In-phase 15-bits codeword
    cq=mod(uq*G,2);              %Quadrature 15-bits codeword
    
    Nbk=(N/k)/2/d; %Number of blocks
    ci_s=[];
    cq_s=[];
    aux=[];
    for p=1:Nbk
        for t=1:n
        aux=ci(((p-1)*d+1):p*d,t)';
        ci_s=[ci_s aux];
        aux=cq(((p-1)*d+1):p*d,t)';
        cq_s=[cq_s aux];
        end
    end
     
    Ac=sqrt(2); %signal amplitude
    c=Ac*((cq_s==0).*(ci_s==0)*(exp(j*(5*pi/4)))+(cq_s==0).*(ci_s==1)...
    *(exp(j*(7*pi/4)))+(cq_s==1).*(ci_s==1)*(exp(j*(9*pi/4)))...
     +(cq_s==1).*(ci_s==0)*(exp(j*(11*pi/4)))); %Transmitted signal with Gray Coding
    
    noise=Ac*sqrt((N0(x)/2))*(randn(1,length(c))+j*randn(1,length(c)));    %AWGN channel
    ray_var=1; %Slow fading coeff. variance
    alpha=sqrt(ray_var*(randn(1,length(c)).^2+(randn(1,length(c)).^2))); %Slow fading coeff. (Rayleigh PDF approximation)

    
    % Variation of fading coeff. every coherence time
    r1 = [ ];
    for p = 1:length(c)/d
        raux = alpha(p)*c(((p-1)*d)+1:(p*d));
        r1 = [r1,raux];
    end
    
    r=r1+noise; %Received signal (Multipath and AWGN over the signal) 
    
    %Seperating bits with same fading
    sr = [];
    for p = 1:length(c)/d
    sraux = r((((p-1)*d)+1:(p*d)))/alpha(p);
    sr = [sr,sraux];
    end
   
    dii=sign(real(sr)); %In-phase hard decision decoding 
    dii(dii<0)=0; %In-phase mapping -1s to 0s again
    dqq=sign(imag(sr)); %Quadrature hard decision decoding
    dqq(dqq<0)=0; %Quadrature mapping -1s to 0s again

    di=[];
    dq=[];
    for p=1:Nbk
        for t=1:n
        di((p-1)*d+1:p*d,t)=dii(1,d*(t-1)+1+6600*(p-1):d*t+6600*(p-1))';
        dq((p-1)*d+1:p*d,t)=dqq(1,d*(t-1)+1+6600*(p-1):d*t+6600*(p-1))';
        end
    end 
  
    syi=mod(di*H',2); %in-phase syndrome calculation
    syq=mod(dq*H',2); %Quadrature syndrome calculation 
    
    e=[zeros(1,n) ; fliplr(eye(n))]; %Error pattern matrix
    
    s_est=mod(e*H',2); %Syndrome matrix to build the look-up table
    
    ci_est=zeros((N/k/2),n); %Initializing in-phase codeword estimation matrix
    cq_est=zeros((N/k/2),n); %Initializing quadrature codeword estimation matrix
    mi=zeros(N/k/2,k); %Initializing in-phase message estimation matrix
    mq=zeros(N/k/2,k); %Initializing quadrature codeword estimation matrix
    
       for p=1:(N/k/2)
           for t=1:(n+1)
               if syi(p,:)==s_est(t,:)
                   ci_est(p,:)=mod(di(p,:)+e(t,:),2); %In-phase estimation of the codeword
                   mi(p,:)=ci_est(p,(n-k+1):n); %In-phase estimation of the message block
               end
               if syq(p,:)==s_est(t,:)
                   cq_est(p,:)=mod(dq(p,:)+e(t,:),2); %Quadrature estimation of the codeword
                   mq(p,:)=cq_est(p,(n-k+1):n); %Quadrature estimation of the message block
               end   
           end
       end
    
    ii=reshape(mi',1,N/2);            %Received in-phase data stream
    qq=reshape(mq',1,N/2);            %Received quadrature data stream
  
    ddata=zeros(1,N); %Decoded data vector initialization
    ddata(1:2:end)=ii; %In-phase hard decision decoding
    ddata(2:2:end)=qq; %Quadrature hard decision decoding
       
    BER(x)=(N-sum(data==ddata))/N; %Calculated BER vector    
end

semilogy(EBN0dB,theory_coded_BER,'g:',EBN0dB,theory_coded_interleaved_BER,'r--',EBN0dB,BER,'o','LineWidth',2)                                
xlabel('E_b/N_0(dB)')                                  
ylabel('BER') 
title('BER of QPSK system with channel coding and interleaving')
legend('Coded-QPSK Theory (Rayleigh and AWGN)',...
    'Coded-QPSK Theory (Rayleigh and AWGN with interleaving)',...
    'Coded-QPSK Simulation (Rayleigh and AWGN with interleaving)');
grid on