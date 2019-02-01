clc
close all;
clear all;                                           

%Simulation of a coded QPSK system in base band simulation
%Channel coding: single error correcting code Hamming (15,11)

n=15; %length of the codeword
k=11; %length of the message
t=1; %Number of corrected bits
N=2285800; %Total number of bits of the data stream 
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

theory_uncoded_BER=(1/2)*erfc(sqrt(EBN0)); %Theoretical BER for uncoded QPSK

Pc=(1/2)*erfc(sqrt(ECN0)); %Theoretical BER for coded QPSK
aux=0; %Auxiliar variable
theory_coded_BER=0;

for p=(t+1):1:n
aux=(1/n)*(p*(factorial(n)/(factorial(p)*factorial(n-p)))*(Pc.^p).*(1-Pc).^(n-p));
theory_coded_BER=theory_coded_BER+aux;
end

for x=1:length(EBN0dB) 
    i=round(rand(1,N));          %In-phase random bit stream
    q=round(rand(1,N));          %Quadrature random bit stream
    
    ui=reshape(i,k,N/k)';        %in-phase 11-bits message
    uq=reshape(q,k,N/k)';        %Quadrature 11-bits message
    
    ci=mod(ui*G,2);              %In-phase 15-bits codeword
    cq=mod(uq*G,2);              %Quadrature 15-bits codeword
    
    ci_s=reshape(ci',1,N*n/k);   %In-phase coded bits stream
    ci_s(ci_s==0)=-1;            %In QPSK 0 = -1
    cq_s=reshape(cq',1,N*n/k);   %Quadrature coded bits stream
    cq_s(cq_s==0)=-1;            %In QPSK 0 = -1
    
    s=ci_s+j*cq_s;   %Transmitted signal in base band
    
    noise=(sqrt(N0(x)/2))*(randn(1,N*n/k)+j*randn(1,N*n/k));    %AWGN channel
    
    r=s+noise;   %Received signal

    di=reshape(sign(real(r)),n,N/k)'; %In-phase hard decision decoding 
    di(di<0)=0;
    dq=reshape(sign(imag(r)),n,N/k)'; %Quadrature hard decision decoding
    dq(dq<0)=0;
    
    syi=mod(di*H',2); %in-phase syndrome calculation
    syq=mod(dq*H',2); %Quadrature syndrome calculation
    
    e=[zeros(1,n) ; fliplr(eye(n))]; %Error pattern matrix
    
    s_est=mod(e*H',2); %Syndrome matrix to build the look-up table
    
    ci_est=zeros((N/k),n); %initializing in-phase codeword estimation matrix
    cq_est=zeros((N/k),n); %initializing quadrature codeword estimation matrix
    mi=zeros(N/k,k); %initializing in-phase message estimation matrix
    mq=zeros(N/k,k); %initializing quadrature codeword estimation matrix
    
       for p=1:(N/k)
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
    
    ii=reshape(mi',1,N);            %Received in-phase data stream
    qq=reshape(mq',1,N);            %Received quadrature data stream

    BER1=(N-sum(i==ii))/N;          %In-phase signal BER calculation
    BER2=(N-sum(q==qq))/N;          %Quadrature signal BER calculation
    BER(x)=mean([BER1 BER2]);       %Total BER     
end

semilogy(EBN0dB,theory_uncoded_BER,'r--',EBN0dB,theory_coded_BER,'b--',EBN0dB, BER,'o-','LineWidth',2)                                
xlabel('E_b/N_0(dB)')                                  
ylabel('BER')
title('BER of QPSK system with (15,11) Hamming code')
legend('Uncoded-QPSK Theory','Coded-QPSK Theory','Coded-QPSK Simulation');
grid on               