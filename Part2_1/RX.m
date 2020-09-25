clear all;
close all;
clc;
tic
%%load('Butter.mat');
%Parameters description
Sampling_Freq = 44100;
% Sampling_Freq = 66207;
Carrier_freq = 1500; %carrier frequency
Nbr_Ofbits=2^10;
Nsc = 128; %number of subcarriers in the OFDM system
Time_of_OFDMs = 58e-3;
Nbr_Of_Samples = round(Time_of_OFDMs*Sampling_Freq); % Number of samples in the OFDM information symbol time
N_Samples = round(Nbr_Of_Samples/Nsc);
disp('record start')
recorder=audiorecorder(Sampling_Freq,16,1);
% Recording time can be changed from here 
recordblocking(recorder,30); 
R=getaudiodata(recorder); 
disp('record stop')
R=R';
figure(1)
plot (R);
t1=0:(Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq)):(length(R)-1)*Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq); 
Carrier = exp(-1i*2*pi*Carrier_freq*t1);
r1= R.*Carrier;
[B,A]=butter(8,0.06);
%Low Pass Filter Filtered_Signal = filter(B,A,r1);
Analog_To_Digital_and_Sampling=Filtered_Signal(1:20:end); %AD conv-downsampling
%%%% Synchronization Finding the normalized periodic pilot(Cauchy-Schwarz formula 
for k=1:length(Analog_To_Digital_and_Sampling)-127
vect1 =Analog_To_Digital_and_Sampling(1*k:k+63);
vect2 =Analog_To_Digital_and_Sampling(k+64:k+127); %
vect2=conj(vect2);
Schwarz(k)=sum(vect1.*vect2);%
(sqrt(sum((abs(vect1).^2)))*sqrt(sum((abs(vect2).^2)))); 
end
figure(3)
plot(abs(Schwarz))
title('Correlated');
%%% Finding the starting point of OFDM signal.
[max_Poit,Pos]=max(abs(Schwarz)); syncronized_rvd_sig_WithPilot=Analog_To_Digital_and_Sampling(Pos-10:end); %(Pos:end)
%%%%%%%%%%%%%%%%%%%% Removing the CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pilot=syncronized_rvd_sig_WithPilot(1:128); % 
Received_Data_With_CP=syncronized_rvd_sig_WithPilot(129:end);
Data_Withot_cp=[];
for k=1:148:length(Received_Data_With_CP)-147
Data_Withot_cp=[Data_Withot_cp Received_Data_With_CP(k+20:k+147)]; %i
end
syncronized_rvd_sig_WithPilot=[pilot Data_Withot_cp];
%%% (Pilot signal generation) FFT and Channel Estimation 
Fast_Forier_Channel=fft(pilot);
x=zeros(1,Nsc);
randn('state', 100);
P = sign(randn(1, Nsc/2)); x(1 : 2 :end) = 2 * P; 

%Channel Estimation for Odd samples
channel_Estimation_uninterpolated=Fast_Forier_Channel./x; % get channel response
%%%%%%%%%%%%%%%%% channel Estimation by Interpolation %%%%%%%%%%%%%%% % we have 1 3 5 ..........127 all odds nbr and we dont have even nbrs but
% we can get it for example by avraging 1 and 3 we get 2 ......
for k=1:2:length(channel_Estimation_uninterpolated)-3
    p=[channel_Estimation_uninterpolated(k) channel_Estimation_uninterpolated(k+2)]; 
    p1=mean(p);
channel_Estimation_uninterpolated(k+1)=p1;
end
channel_est_interpolated=channel_Estimation_uninterpolated; %now we have the signal data+channel response
%%%%%%%% We are removing the effect of the channel %%%%%%
Pure_Data=[];
for k=1:128:length(syncronized_rvd_sig_WithPilot)-127
QQ=fft(syncronized_rvd_sig_WithPilot(k:k+127),128);
Pure_Data=[Pure_Data QQ./channel_est_interpolated]; 
end
%now we got the actual signal data 
Pure_Data_Witho_Pilot = Pure_Data(129:end);
%%%%%%%%%%%%%%% DEMAPPING %%%%%%%%%%%%%%%%%%%%%%%%%%%
Received_bits=[];
for k=1:1:length(Pure_Data_Witho_Pilot)
Real = real(Pure_Data_Witho_Pilot); 
Img = imag(Pure_Data_Witho_Pilot); 
if Real(k)>0 && Img (k)>0
Received_bits(2*k-1)=1;
Received_bits(2*k)=1; 
elseif Real(k)>0 && Img (k)< 0
Received_bits(2*k-1)=1;
Received_bits(2*k)=0;
elseif Real(k)< 0 && Img (k)>0 
    Received_bits(2*k-1)=0;
    Received_bits(2*k)=1;
elseif Real(k)<=0 && Img (k)<=0
    Received_bits(2*k-1)=0; 
    Received_bits(2*k)=0;
end
end
%constLen=6; % specifying the delay for each of the k input bit streams. %trellis=poly2trellis(6,[77,45]);
% message_bits=vitdec( Received_bits,trellis,50,'term','hard');
%
%Viterbi decoding of Data
%1.ok
% trellis=poly2trellis(6,[77,45]);
% tblen=50;
% Uncoded_receved_Data=vitdec(Received_bits,trellis,50,'term','hard'); & % 1.ok
constlen=6;
%codegen = [77 45]; % Polynomial
trellis = poly2trellis(constlen, [77 45]);
Uncoded_receved_Data=vitdec(Received_bits,trellis,50,'term','hard');
%Uncoded_receved_Data=Uncoded_receved_Data(1:end); % We will get the massege if it increas massage + nonces text
% Finding the length of the original masseage %%%1.ok FOR TEXT
% Msg_L=Uncoded_receved_Data(1:15); %
%
% Msg_Len=[];
% for i=1:length(Msg_L)
% Msg_Len=[Msg_Len mat2str(Msg_L(i))];
% end %
% Msg_Len_nbr=bin2dec(Msg_Len);
% Uncoded_Data=Uncoded_receved_Data(129:129+Msg_Len_nbr-1);
% % % [n1,r1] = biterr(Received_bits(tblen+1:end),Uncoded_Data(1:end-tblen));
% % % disp([' Error rates are: ',num2str(r1)]) %
% wh=2.^[6:-1:0];
% Secret_Massege=char(Uncoded_Data(1:7)*wh');
% for l=2:floor(length(Uncoded_Data)/7)
% Secret_Massege=[Secret_Massege char(Uncoded_Data(7*(l-1)+1:7*l)*wh')];
% end
% Secret_Massege %%% 1.ok
%%%%2 ok FOR IMAGE
Massege_Len=Uncoded_receved_Data(1:30); Massg_Len=[];
for i=1:length(Massege_Len)
Massg_Len=[Massg_Len mat2str(Massege_Len(i))]; 
end

Msg_Len_nbr=bin2dec(Massg_Len);
%%%% Finding nbr of Rows Massege_Len2=Uncoded_receved_Data(31:60); Row_len=[];
for i=1:length(Massege_Len2)
Row_len=[Row_len mat2str(Massege_Len2(i))];
end
Nbr_Rows=bin2dec(Row_len);
Colu_Len=Uncoded_receved_Data(61:90); Colum_len=[];
for i=1:length(Colu_Len)
Colum_len=[Colum_len mat2str(Colu_Len(i))];
end
Nbr_Colums=bin2dec(Colum_len);
Uncoded_Data=Uncoded_receved_Data(129:129+Msg_Len_nbr-1);
Convert_Matrix= reshape(Uncoded_Data,[Nbr_Rows,Nbr_Colums]); figure(1)
imshow(Convert_Matrix)


