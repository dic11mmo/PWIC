%%TASK 2
clc
clear all
close all
%Parameters description
Sampling_Freq = 44100; Carrier_Freq = 4000;
Ts = 2.3*10^-3;
Beta = (1/100);
sample_time = 1/Sampling_Freq;
N = round(Ts*Sampling_Freq); Nbr_Ofbits=1000;
mu=0.00023;
%Generate half cycle sinus Pulse shape:
t=0:sample_time:((N-1)/Sampling_Freq); puls_shap=sin((pi/Ts)*t);
%Bits generation
Transmit_bits=(rand(1,Nbr_Ofbits)>0.5); vec_1Andminusone=2*Transmit_bits-1;;
% y is a vector of 1 and -1
Real=vec_1Andminusone(1:2:Nbr_Ofbits-1); Img=vec_1Andminusone(2:2:Nbr_Ofbits);
Imag_part =1i*Img ;
Complex_sig=( Real + Imag_part); % QAM Symbols
%Generate carrier frequency and pulse train
pulse_train=[]; %Train pulse (empty vector) 
for k=1:1:length(Complex_sig);
pulse_train=[pulse_train (puls_shap*Complex_sig(k))];
end
% Real and img part of the Pulse train
Real_PulsTrain= real(pulse_train); Img_PulsTrain = imag(pulse_train);
%%%%%%%%%%                      Modulation     %%%%%
t1=0:1/Sampling_Freq: ((length(Real_PulsTrain)-1)/Sampling_Freq); RQ=sqrt(2).*Real_PulsTrain.*cos(2*pi*Carrier_Freq*t1); RI=-Img_PulsTrain.*sqrt(2).*sin(2*pi*Carrier_Freq*t1);
Modulat= RQ + RI;
% To find nbr of delays
Differe_TwoSignas = round(Sampling_Freq*mu);
for i =Differe_TwoSignas +1:(N* Nbr_Ofbits/2)
Matrix_r(i)=sqrt((1-Beta^2))*Modulat(i) + Beta*Modulat(i- Differe_TwoSignas);
end
figure(4)
plot(Matrix_r)
M = 2;
Es=sum((abs(Complex_sig).^2)/Sampling_Freq); Eb=Es/log2(M); % Eb =Es/k 2.31 digital communication 
SNR_dB= 1:.5:15;
BER = zeros(1,length(SNR_dB));
SNR_Linear = (10.^(SNR_dB/10));
x = 1;
for l=(SNR_Linear);
N0=Eb./(l);
%Definition of noise
sigma=sqrt(Sampling_Freq*(N0/2));
Nois=sigma.*randn(1,length(Matrix_r));
%Channel with noise
Signal_With_Noise=Matrix_r+Nois;
%%%%%   QPSK_Demodulation %%%%%%%%%
RQ = Signal_With_Noise.*cos(2*pi*Carrier_Freq*t1); RI =-Signal_With_Noise.*(sin(2*pi*Carrier_Freq*t1)); Demodulat = RQ +RI;
%%%%%%%%%% https://www.youtube.com/watch?v=h46MLgZDQqU %%%%%%%%%% % Matching filter
%Matched filter is the convolution between r and Pulse shape:
Match_filter_RQ=conv(RQ,puls_shap); 
Match_filter_RI=conv(RI,puls_shap);
%Sampling the received signal and detection
% we get top of the signal for sampling 
Rec_signal_Qpart=Match_filter_RQ(1:N:end);
% we get a vector of values from 1 then every 101:e
% sign gives us +-1 vectors we start from 2
% because we have zero in the starting of the vector 
Rec_signal_Qpart_TopPoints = sign(Rec_signal_Qpart(2:1:end));
Rec_signal_Ipart=(Match_filter_RI(1:N:end)); Rec_signal_Ipart_Img=1i*Rec_signal_Ipart;
% sign gives us +-1 vectors we start from 2
% because we have zero in the starting of the vector
Rec_signal_Ipart_TopPoints = sign(Rec_signal_Ipart_Img(2:1:end));
Received_Complex = Rec_signal_Qpart_TopPoints + Rec_signal_Ipart_TopPoints ;
% b=1
Received_bits=zeros(1,100);
for k=1:length(Received_Complex)
if Received_Complex(k)==1+1i; 
    %Real(k) >= 0 && Img(k) >= 0
    Received_bits(2*k-1)=0; 
    Received_bits(2*k)=0;
elseif Received_Complex(k)==-1+1i;
    % Real(k) < 0 && Img(k) > 0
Received_bits(2*k-1)=1;
Received_bits(2*k)=0;
elseif Received_Complex(k)==-1-1i;
        %Real(k) < 0 && Img(k) < 0


Received_bits(2*k-1)=1;
Received_bits(2*k)=1;
elseif Received_Complex(k)==1-1i;
        %Real(k) > 0 && Img(k) < 0
Received_bits(2*k-1)=0; Received_bits(2*k)=1;
end
end
Received_bits2 = [];
for j = 1 : length(Received_bits)
if Received_bits(j)==1 ;
Received_bits2=[Received_bits2 (Received_bits(j)*0)];
elseif Received_bits(j)==0;
Received_bits2=[Received_bits2 (Received_bits(j)+1)];
end
end
%BER = (nbr of errors)/(Total nbr of bits sent)
% BER(x)=length(find(Received_bits2-Transmit_bits))/Nbr_Ofbits;
BER(x)=length(find(Received_bits2-Transmit_bits))/Nbr_Ofbits; x=x+1;
end
tBer_QAM = (1/4)*3*1/2*erfc(sqrt(4*0.05*(10.^(SNR_dB/10)))); % theoretical ber of QAM
figure(1);
plot(t,puls_shap);% Generate Sine Wave
title(' Half cycle sinus Base pulse')
figure(2)
subplot(2,1,1)
plot(Real_PulsTrain, 'r')
title('Real part of signal altenatives s(t) ') 
ylabel('Amplitude')
xlabel( 'Signal time t ')
hold on
subplot(2,1,2)
plot(Img_PulsTrain,'g')
%stem(Img_PulsTrain)
title('Imaginare part of signal altenatives s(t)');
figure(3);
 subplot(2,1,1)
plot(Modulat);
title(' Modulated signal')
subplot(2,1,2)
plot(Demodulat);
title(' Demodulated signal')
figure(4);
subplot(2,1,1)
plot( Matrix_r);
title(' Channel without noise')
hold on
subplot(2,1,2)
plot( Signal_With_Noise);
title(' Channel with white Gaussian Noise') % figure(5)
% semilogy(SNR_dB ,BER);
% title('BER function of different SNR');
% xlabel('Eb/N0');
% ylabel('BER');
figure(5)
semilogy(SNR_dB,tBer_QAM,'-*r','linewidth',3); hold on
semilogy(SNR_dB ,BER);
grid on
legend('Theoretical QAM',' Our result');
 title('BER function of different SNR');
xlabel('Eb/N0');
ylabel('BER');


