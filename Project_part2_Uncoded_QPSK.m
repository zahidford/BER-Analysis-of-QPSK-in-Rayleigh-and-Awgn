clc
close all;
clear all;

%Simulation of uncoded QPSK system in base band (AWGN + Rayleigh fading channel)

N=1e6; %Total number of bits of the data stream
EBN0dB=0:1:20; %Normalized SNR per bit in dB
EBN0=10.^(EBN0dB/10); %Normalized SNR per bit
N0=10.^(-EBN0dB/10); %Noise spectral density

theory_BER_AWGN=(1/2)*erfc(sqrt(EBN0)); %Theoretical uncoded QPSK BER (AWGN)
theory_BER_Rayleigh=(1/2)*(1-sqrt(EBN0./(1+EBN0))); %Theoretical uncoded QPSK BER with Rayleign fading

for x=1:length(EBN0dB) 
    data=round(rand(1,N)); %Random data stream

    i=data(1:2:end); %In-phase bits
    q=data(2:2:end); %Quadrature bits

    Ac=sqrt(2); %signal amplitude
    s=Ac*((q==0).*(i==0)*(exp(j*(5*pi/4)))+(q==0).*(i==1)...
    *(exp(j*(7*pi/4)))+(q==1).*(i==1)*(exp(j*(9*pi/4)))...
     +(q==1).*(i==0)*(exp(j*(11*pi/4)))); %Transmitted signal with Gray Coding
 
    noise=Ac*sqrt(N0(x)/2)*(randn(1,N/2)+j*randn(1,N/2)); %AWGN channel
    ray_var=1; %slow fading coeff. variance
    alpha=sqrt(ray_var*((randn(1,N/2)).^2+(randn(1,N/2)).^2)); %Slow fading coeff. (Rayleigh PDF approximation)

    r=alpha.*s+noise; %Received signal (Multipath and AWGN over the signal)

    ddata=zeros(1,N); %Decoded data vector initialization
    ddata(1:2:end)=sign(real(r)); %In-phase hard decision decoding
    ddata(2:2:end)=sign(imag(r)); %Quadrature hard decision decoding
    
    ddata(ddata==-1)=0; %Mapping -1s to 0s again
    
    BER(x)=(N-sum(data==ddata))/N; %Calculated BER vector
end

semilogy(EBN0dB,theory_BER_AWGN,'g:',EBN0dB,theory_BER_Rayleigh,'r--',EBN0dB, BER,'o','LineWidth',2)                                
xlabel('E_b/N_0(dB)')                                  
ylabel('BER') 
title('BER of QPSK system without channel coding')
legend('Uncoded-QPSK Theory (AWGN)',...
    'Uncoded-QPSK Theory (Rayleigh and AWGN)',...
    'Uncoded-QPSK Simulation (Rayleigh and AWGN)');
grid on