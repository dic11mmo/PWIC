%%TASK 3
clc
clear all
close all
%Parameters description
load('signal3.mat');
%%load('Butter.mat');
%%%%%%%%%%%%%%%%%%%%%%% SYSTEM PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Sampling_Freq=44100;
Tsamp=1/Sampling_Freq;
Carrier_freq=10000; % carrier frequency
Nbr_Of_sub=128; %number of subcarriers in the OFDM system
Length_cyclic_prefix=20; %length of cyclic prefix 
Time_of_OFDMs=58*10^(-3); 
Time_subcarriers=Time_of_OFDMs/Nbr_Of_sub;
Samp_freq_Of_Subcarrie=1/Time_subcarriers; 
% Subcarrier separation c =round(Sampling_Freq/Samp_freq_Of_Subcarrie); OFDM_sample_rate=Time_subcarriers/Tsamp;
Duration_Cyclic_prefix = (58*10^-3)*(Length_cyclic_prefix/Nbr_Of_sub); TsCp = Time_of_OFDMs + Duration_Cyclic_prefix;
Carrier = exp(-j*2*pi*Carrier_freq*t); r1=R.*Carrier(1:length(R));
[B,A]=butter(8,0.05);
%Low Pass Filter
Filtered_Signal = filter(B,A,r1); Analog_To_Digital_and_Sampling=Filtered_Signal(1:20:end); %AD conv- downsampling
%%%% Synchronization Finding the normalized periodic pilot(Cauchy-Schwarz formula
for k=1:length(Analog_To_Digital_and_Sampling)-127 vect1 =Analog_To_Digital_and_Sampling(1*k:k+63); vect2 =Analog_To_Digital_and_Sampling(k+64:k+127);
% Schwarz(i)=sum(vect1.*conj(vect2)/(sqrt(sum(abs(vect1).^2)))*sqrt(sum(abs(con j(vect2).^2))));
vect2=conj(vect2);
Schwarz(k)=sum(vect1.*vect2)/(sqrt(sum((abs(vect1).^2)))*sqrt(sum((abs(vect2) .^2))));
end
figure(2)
plot(abs(Schwarz))
title('Correlated');
%%% Finding the starting point of OFDM signal.
[max Pos]=max(abs(Schwarz)); 
syncronized_received_signal=Analog_To_Digital_and_Sampling(Pos-10:end);
%%%%%%%%%%%%%%%%%%%% Removing the CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pilot=syncronized_received_signal(1:128); %
Received_Data_With_CP=syncronized_received_signal(129:end);
Remove_Cp=[];
for k=1:148:length(Received_Data_With_CP)-147
Remove_Cp=[Remove_Cp Received_Data_With_CP(k+20:k+147)];
%i
end
syncronized_received_signal=[pilot Remove_Cp];
%%% (Pilot signal generation) FFT and Channel Estimation
Fast_Forier_Channel=fft(pilot); x=zeros(1,Nbr_Of_sub); randn('state', 100);
P = sign(randn(1, Nbr_Of_sub/2)); x(1 : 2 :end) = 2 * P;
%Channel Estimation for Odd samples
channel_Estimation_uninterpolated=Fast_Forier_Channel./x; % get channel response
%%%%%%%%%%%%%%%%% channel Estimation by Interpolation %%%%%%%%%%%%%%%
% we have 1 3 5 ..........127 all odds nbr and we dont have even nbrs but % we can get it for example by avraging 1 and 3 we get 2 ......
for k=1:2:length(channel_Estimation_uninterpolated)-2 
    p=[channel_Estimation_uninterpolated(k);
channel_Estimation_uninterpolated(k+2)]; p1=mean(p);
channel_Estimation_uninterpolated(k+1)=p1;
end
channel_est_interpolated=channel_Estimation_uninterpolated;
%now we have the signal data+channel response
%%%%%%%% We are removing the effect of the channel %%%%%%
Pure_Data=[];
for k=1:128:length(syncronized_received_signal)-127
QQ=fft(syncronized_received_signal(k*1:k+127),128);


Pure_Data=[Pure_Data QQ./channel_est_interpolated];
end
%now we got the actual signal data
%%%%%%%%%%%%%%% DEMAPPING %%%%%%%%%%%%%%%%%%%%%%%%%%%
Received_bits=[];
for k=1:1:length(Pure_Data)
Real = real(Pure_Data); Img = imag(Pure_Data);
if Real(k)>0 && Img (k)>0
Received_bits(2*k-1)=0;
Received_bits(2*k)=0; elseif Real(k)>0 && Img (k)< 0
Received_bits(2*k-1)=0;
Received_bits(2*k)=1; elseif Real(k)< 0 && Img (k)>0
Received_bits(2*k-1)=1;
Received_bits(2*k)=0; elseif Real(k)<0 && Img (k)<0
Received_bits(2*k-1)=1;
Received_bits(2*k)=1;
end
end
Rec_bits_removepilot=Received_bits(257:end);
length_message_bits=Rec_bits_removepilot(1:256);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Getting the length of the message %%%%%%%%%%%%%%%%%%%%%%
% poly2trellis : converts a polynomial representation of a feedforward %convolutional encoder to a trellis structure. For a rate k/n code,
% the encoder input is a vector of and length k, and the encoder
% output is a vector of length n. Therefore,
constLen=6; % specifying the delay for each of the k input bit streams. 
trellis=poly2trellis(6,[77,45]);
message_bits=vitdec(length_message_bits,trellis,50,'term','hard');
%bi2de Convert binary vectors to decimal numbers
Length_Of_Recev_Massege=bi2de(message_bits(1:10));
Coded_receved_data=Rec_bits_removepilot(257:end);
%Viterbi decoding of Data
Uncoded_receved_Data=vitdec(Coded_receved_data,trellis,30,'term','hard');
% We will get the massege if it increas massage + nonces text by changing Uncoded_receved_Data(1:602) increase 
% or decrease length of data 602 up or down you will get end of the data.
Uncoded_receved_Data=Uncoded_receved_Data(1:602); 
% ASCII
wh=2.^[6:-1:0]; Secret_Massege=char(Uncoded_receved_Data(1:7)*wh');
for l=2:floor(length(Uncoded_receved_Data)/7), 
    Secret_Massege=[Secret_Massege char(Uncoded_receved_Data(7*(l- 1)+1:7*l)*wh')];
end
Secret_Massege
