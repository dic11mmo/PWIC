

function [CRC,Received_packet, Accept_Or_Reject ] = RX( recorded_signal )
%  [ Received_packet, Accept_Or_Reject ] = RX( recorded_signal )
%[Received_packet,Accept_Or_Reject ]
load('Butter.mat');
 load('randomvector1.mat');
%Parameters description
Sampling_Freq = 44100;
% Sampling_Freq = 66207;
Carrier_freq = 3000;  %carrier frequency

Nsc = 128; %number of subcarriers in the OFDM system
Time_of_OFDMs = 58e-3;
Nbr_Of_Samples = round(Time_of_OFDMs*Sampling_Freq); % Number of samples in the OFDM information symbol time
N_Samples = round(Nbr_Of_Samples/Nsc);
 
 R = recorded_signal'; 
 t1=0:(Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq)):(length(R)-1)*Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq);
Carrier = exp(-1i*2*pi*Carrier_freq*t1);
 r1= R.*Carrier;
 

[B,A]=butter(8,0.06);
%Low Pass Filter
Filtered_Signal = filter(B,A,r1);
 
Analog_To_Digital_and_Sampling=Filtered_Signal(1:20:end); %AD conv-downsampling

%%%% Synchronization Finding the normalized periodic pilot(Cauchy-Schwarz formula
for k=1:length(Analog_To_Digital_and_Sampling)-127  
    vect1 =Analog_To_Digital_and_Sampling(1*k:k+63);
    vect2 =Analog_To_Digital_and_Sampling(k+64:k+127);
% 
    vect2=conj(vect2);
    Schwarz(k)=sum(vect1.*vect2);%(sqrt(sum((abs(vect1).^2)))*sqrt(sum((abs(vect2).^2))));
end
figure(3)
plot(abs(Schwarz))
title('Correlated');
%%% Finding the starting point of OFDM signal.
 [max_Poit,Pos]=max(abs(Schwarz));
syncronized_rvd_sig_WithPilot=Analog_To_Digital_and_Sampling(Pos-10:end); %(Pos:end)
 
%%%%%%%%%%%%%%%%%%%%    Removing the CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pilot=syncronized_rvd_sig_WithPilot(1:128); % 
Received_Data_With_CP=syncronized_rvd_sig_WithPilot(129:end);

Data_Withot_cp=[];
for k=1:148:length(Received_Data_With_CP)-147
   Data_Withot_cp=[Data_Withot_cp Received_Data_With_CP(k+20:k+147)]; 
%    i
end
syncronized_rvd_sig_WithPilot=[pilot Data_Withot_cp];

%%% (Pilot signal generation) FFT and Channel Estimation
Fast_Forier_Channel=fft(pilot);
x=zeros(1,Nsc);
randn('state', 100);
P = sign(randn(1, Nsc/2)); 
x(1 : 2 :end) = 2 * P; 
%Channel Estimation for Odd samples
channel_Estimation_uninterpolated=Fast_Forier_Channel./x; % get channel response

%%%%%%%%%%%%%%%%% channel Estimation by Interpolation %%%%%%%%%%%%%%%
% we have 1 3 5 ..........127 all odds nbr and we dont have even nbrs but
% we can get it for example by avraging 1 and 3 we get 2 ......
    for k=1:2:length(channel_Estimation_uninterpolated)-3
        p=[channel_Estimation_uninterpolated(k) channel_Estimation_uninterpolated(k+2)];
        p1=mean(p); 
        channel_Estimation_uninterpolated(k+1)=p1;
    end
channel_est_interpolated=channel_Estimation_uninterpolated;
%now we have the signal data+channel response
 
%%%%%%%% We are removing the effect of the channel  %%%%%%
Pure_Data=[]; 
for k=1:128:length(syncronized_rvd_sig_WithPilot)-127
    QQ=fft(syncronized_rvd_sig_WithPilot(k:k+127),128);
    Pure_Data=[Pure_Data QQ./channel_est_interpolated];
end
%now we got the actual signal data
Pure_Data_Witho_Pilot = Pure_Data(129:end);

%%%%%%%%%%%%%%%     DEMAPPING      %%%%%%%%%%%%%%%%%%%%%%%%%%%
Received_bits=[];
for k=1:1:length(Pure_Data_Witho_Pilot) 
    Real = real(Pure_Data_Witho_Pilot);
    Img = imag(Pure_Data_Witho_Pilot);
    if  Real(k)>0 && Img (k)>0
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
   

trellis = poly2trellis(6,[77,45]);   % constraint length is 6
Un_Coded_Data = vitdec(Received_bits,trellis,50,'trunc','hard');
Un_CodeData_With_CRC= Un_Coded_Data(1:1024);
    Un_CodeData_WithOutCRC1 = xor(randomvector1,Un_CodeData_With_CRC);

Received_packet =Un_CodeData_WithOutCRC1(1:end-16); 


%%%CRC detector

    detect = comm.CRCDetector([16 15 2 0],'ChecksumsPerFrame',1);
   
   [~, CRC] = step(detect,Un_CodeData_WithOutCRC1');
    Accept_Or_Reject=[];
   Random_Nbr = randi([0,1],1,896);
    % If CRC = 1 than there is erros and Rejected 
    if (CRC == 0) 
        Accept_Or_Reject=zeros(1,128);
        Accept_Or_Reject = [Accept_Or_Reject Random_Nbr];
        Ack_or_Nack=true;
        %If CRC = 0 than there is no erros and accept it 
    elseif (CRC == 1)
        Accept_Or_Reject=ones(1,128);
        Accept_Or_Reject = [Accept_Or_Reject Random_Nbr];
        Ack_or_Nack=false;
    end







%     save Accept_Or_Reject
%     
%     %Parameters description
%     load('Accept_Or_Reject');
% % Sampling_Freq = 44100;
% % Carrier_freq = 000;    %carrier frequency
% % %Nbr_Ofbits=2^10;
% % Nsc = 128; %number of subcarriers in the OFDM system
% % Time_of_OFDMs = 58e-3;
% % Nbr_Of_Samples = round(Time_of_OFDMs*Sampling_Freq); % Number of samples in the OFDM information symbol time
% % N_Samples = round(Nbr_Of_Samples/Nsc);
%  
% 
%  % to make length of data multiply of 128 
% 
% % gen = comm.CRCGenerator([16 15 2 0],'ChecksumsPerFrame',1);
% % codeword = step(gen, Accept_Or_Reject');
% % codeword =codeword';
%   codeword = Accept_Or_Reject ;
% %bits=(rand(1,Nbr_Ofbits)>0.5);
% %bits=rand(1,Nbr_Ofbits)>0.5;
% %Apply convolutionnal code
%  trellis = poly2trellis(6,[77 , 45]);
%  bits = convenc(codeword,trellis);
%  
% 
%  
% 
% %Mapping between input_01 and symbol 
%  
% 
% vec_1Andminusone=2*bits-1; % y is a vector of 1 and -1
%  
% 
%  
% 
% Re=vec_1Andminusone(1:2:length(bits)-1); 
% Im=vec_1Andminusone(2:2:length(bits));
% Complex_sig=Re+1i*Im; % QAM Symbols 
% 
% Adding_Pilots = adding_pilots(Nsc)';
% % 3. Adding the pilot at the beginning of the symbols
% symbols=[Adding_Pilots Complex_sig];
% % 4.Converting Serial-to-parallel by using reshape commond
% Serial_to_parallel=reshape(symbols,[Nsc,length(symbols)/Nsc]);
% % 5. Taking the IFFT of the 4
% Result_Of_ifft=ifft(Serial_to_parallel,Nsc);
% % 6.Converting Parallel to Serial by using reshape commond
% Paralle_To_Serial=reshape(Result_Of_ifft,[1,Nsc*length(symbols)/Nsc]);
%  %Remove_zeros = Removing_Paddad_Zeros(Zero_Paddad_Data)
% %7.Finding the CPs which is the last 20 values of each Nsc placed at the
% %beginning)
% Pilot = Paralle_To_Serial(1:1:128);
% Paralle_To_Serial_WithPilot= Paralle_To_Serial(129:1:end);
% for l=1:length(Paralle_To_Serial_WithPilot)/Nsc
%     CP(l,:)=Paralle_To_Serial_WithPilot(l*Nsc-19:1:l*Nsc);
% end
% %8.The CPs inserts
% Adding_Cp_WithData=Paralle_To_Serial_WithPilot;
% Adding_Cp_WithData=[CP(1,:) Adding_Cp_WithData];
% for m=1:1:length(CP(:,1))-1
%     Adding_Cp_WithData=[Adding_Cp_WithData(1:1:m*Nsc+m*20) CP(m+1,:) Adding_Cp_WithData(m*Nsc+m*20+1:1:end)];
% end 
% %9.Adding pilot again
%  Add_Pilot_Data_WithCp =[ Pilot Adding_Cp_WithData] ; 
% %10.Digital to Analog convertion 
% %Digital_To_Analog=interp1(1:1:length(Add_Pilot_Data_WithCp), Add_Pilot_Data_WithCp, 1:1/22:length(Adding_Cp_WithData));
% Digital_To_Analog=interp(Add_Pilot_Data_WithCp,20);
% %Digital_To_Analog=Digital_To_Analog(1:end-19);
% % figure(1)
% % plot(abs(Digital_To_Analog))
% % Modulation 
% Real_Part_Ofsig = real(Digital_To_Analog);
% Imag_Part_Ofsig = imag(Digital_To_Analog);
%  
% 
%  
% 
% %t1=0:(1/Sampling_Freq):(length(Real_Part_Ofsig)-1)/Sampling_Freq;
% t1=0:(Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq)):(length(Real_Part_Ofsig)-1)*Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq);
% Carrier = exp(1i*2*pi*Carrier_freq*t1);
% Transmitting_Signa=Digital_To_Analog.*Carrier;
%  
%  %Transmitting_Signa=[ Transmitting_Signa]; 
% Transmitting_Signa=real( Transmitting_Signa);
%  figure(2)
%  plot(Transmitting_Signa)
% title(' Trasnmitted signal');
% 
% %Transmitting via microphone 
%  Sending=audioplayer(Transmitting_Signa,Sampling_Freq); 
% playblocking(Sending);


%   

























 




















% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % clear;
% % close all;
% % clc;
% %  
% % 
% % function [Received_Bits, Accept_Or_Reject] = RX_Transmitt( recorded_signal )
% % load('Butter.mat');
% %Parameters description
% Sampling_Freq = 44100;
% % Sampling_Freq = 66207;
% Carrier_freq = 1500;  %carrier frequency
% Nbr_Ofbits=2^10;
% Nsc = 128; %number of subcarriers in the OFDM system
% Time_of_OFDMs = 58e-3;
% Nbr_Of_Samples = round(Time_of_OFDMs*Sampling_Freq); % Number of samples in the OFDM information symbol time
% N_Samples = round(Nbr_Of_Samples/Nsc);
%  
% 
%  
% 
% disp('record start')
% recorder=audiorecorder(Sampling_Freq,16,1);
% recordblocking(recorder,30);
%  
% 
% R=getaudiodata(recorder); 
% disp('record stop')
%  
% 
% R=R';
% figure(1)
% plot (R);
%  t1=0:(Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq)):(length(R)-1)*Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq);
% Carrier = exp(-1i*2*pi*Carrier_freq*t1);
%  r1= R.*Carrier;
%  
% 
% [B,A]=butter(8,0.06);
% %Low Pass Filter
% Filtered_Signal = filter(B,A,r1);
%  
% 
%  
% 
% Analog_To_Digital_and_Sampling=Filtered_Signal(1:20:end); %AD conv-downsampling
%  
% 
%  
% 
% %%%% Synchronization Finding the normalized periodic pilot(Cauchy-Schwarz formula
% for k=1:length(Analog_To_Digital_and_Sampling)-127  
%     vect1 =Analog_To_Digital_and_Sampling(1*k:k+63);
%     vect2 =Analog_To_Digital_and_Sampling(k+64:k+127);
% % 
%     vect2=conj(vect2);
%     Schwarz(k)=sum(vect1.*vect2);%(sqrt(sum((abs(vect1).^2)))*sqrt(sum((abs(vect2).^2))));
% end
% figure(3)
% plot(abs(Schwarz))
% title('Correlated');
% %%% Finding the starting point of OFDM signal.
%  [max_Poit,Pos]=max(abs(Schwarz));
% syncronized_rvd_sig_WithPilot=Analog_To_Digital_and_Sampling(Pos-10:end); %(Pos:end)
%  
% 
%  
% 
% %%%%%%%%%%%%%%%%%%%%    Removing the CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pilot=syncronized_rvd_sig_WithPilot(1:128); % 
% Received_Data_With_CP=syncronized_rvd_sig_WithPilot(129:end);
%  
% 
%  
% 
% Data_Withot_cp=[];
% for k=1:148:length(Received_Data_With_CP)-147
%    Data_Withot_cp=[Data_Withot_cp Received_Data_With_CP(k+20:k+147)]; 
% %    i
% end
% syncronized_rvd_sig_WithPilot=[pilot Data_Withot_cp];
%  
% 
%  
% 
% %%% (Pilot signal generation) FFT and Channel Estimation
% Fast_Forier_Channel=fft(pilot);
% x=zeros(1,Nsc);
% randn('state', 100);
% P = sign(randn(1, Nsc/2)); 
% x(1 : 2 :end) = 2 * P; 
% %Channel Estimation for Odd samples
% channel_Estimation_uninterpolated=Fast_Forier_Channel./x; % get channel response
%  
% 
%  
% 
% %%%%%%%%%%%%%%%%% channel Estimation by Interpolation %%%%%%%%%%%%%%%
% % we have 1 3 5 ..........127 all odds nbr and we dont have even nbrs but
% % we can get it for example by avraging 1 and 3 we get 2 ......
%     for k=1:2:length(channel_Estimation_uninterpolated)-3
%         p=[channel_Estimation_uninterpolated(k) channel_Estimation_uninterpolated(k+2)];
%         p1=mean(p); 
%         channel_Estimation_uninterpolated(k+1)=p1;
%     end
% channel_est_interpolated=channel_Estimation_uninterpolated;
% %now we have the signal data+channel response
%  
% 
%  
% 
% %%%%%%%% We are removing the effect of the channel  %%%%%%
%  
% 
%  
% 
% Pure_Data=[]; 
% for k=1:128:length(syncronized_rvd_sig_WithPilot)-127
%     QQ=fft(syncronized_rvd_sig_WithPilot(k:k+127),128);
%     Pure_Data=[Pure_Data QQ./channel_est_interpolated];
% end
% %now we got the actual signal data
% Pure_Data_Witho_Pilot = Pure_Data(129:end);
%  
% 
% %%%%%%%%%%%%%%%     DEMAPPING      %%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
% 
%  
% 
% Received_bits=[];
% for k=1:1:length(Pure_Data_Witho_Pilot) 
%     Real = real(Pure_Data_Witho_Pilot);
%     Img = imag(Pure_Data_Witho_Pilot);
%     if  Real(k)>0 && Img (k)>0
%     
% 
%  
% 
%         Received_bits(2*k-1)=1;
%         Received_bits(2*k)=1;
%     elseif Real(k)>0 && Img (k)< 0
%         Received_bits(2*k-1)=1;
%         Received_bits(2*k)=0;
%     elseif Real(k)< 0 && Img (k)>0
%         Received_bits(2*k-1)=0;
%         Received_bits(2*k)=1;
%     elseif Real(k)<=0 && Img (k)<=0
%         Received_bits(2*k-1)=0;
%         Received_bits(2*k)=0;
%     end
% end
%  
% 
%  %CRC detector
%     detect = comm.CRCDetector([1 0 0 1],'ChecksumsPerFrame',1);
%     [~, CRC] = step(detect,Rbits')
%     
% 
%     
% 
%     
% 
%     
% 
%     Accept_Or_Reject=[];
%     % If CRC = 1 than there is erros and Rejected 
%     if (CRC == 1)
%         Accept_Or_Reject=ones(1,253);
%         Ack_or_Nack=false;
%         If CRC = 0 than there is no erros and accept it 
%     elseif (CRC == 0) 
%         Accept_Or_Reject=zeros(1,(2*Nsc)-3);
%         Ack_or_Nack=true;
%     end
%  
% 
%     save Accept_Or_Reject
%     
% 
% %%%%2 ok FOR IMAGE 
% % 
% % Massege_Len=Uncoded_receved_Data(1:30);
% % Massg_Len=[];
% % for i=1:length(Massege_Len)
% %     Massg_Len=[Massg_Len mat2str(Massege_Len(i))];
% % end
% % 
% % Msg_Len_nbr=bin2dec(Massg_Len);
% % 
% % %%%% Finding nbr of Rows
% % Massege_Len2=Uncoded_receved_Data(31:60);
% % Row_len=[];
% % for i=1:length(Massege_Len2)
% %     Row_len=[Row_len mat2str(Massege_Len2(i))];
% % end
% % 
% % Nbr_Rows=bin2dec(Row_len);
% % 
% % Colu_Len=Uncoded_receved_Data(61:90);
% % Colum_len=[];
% % for i=1:length(Colu_Len)
% %     Colum_len=[Colum_len mat2str(Colu_Len(i))];
% % end
% % 
% % Nbr_Colums=bin2dec(Colum_len);
% % 
% % %  
% % Uncoded_Data=Uncoded_receved_Data(129:129+Msg_Len_nbr-1);
% %  Convert_Matrix= reshape(Uncoded_Data,[Nbr_Rows,Nbr_Colums]);
% %  figure(1)
% %   imshow(Convert_Matrix)
%  
% 
%  
% 
% %Parameters description
% Sampling_Freq = 44100;
% Carrier_freq = 1500;    %carrier frequency
% %Nbr_Ofbits=2^10;
% Nsc = 128; %number of subcarriers in the OFDM system
% Time_of_OFDMs = 58e-3;
% Nbr_Of_Samples = round(Time_of_OFDMs*Sampling_Freq); % Number of samples in the OFDM information symbol time
% N_Samples = round(Nbr_Of_Samples/Nsc);
%  
% 
% load (Accept_Or_Reject)
%  %CRC check 
%     gen = comm.CRCGenerator([1 0 0 1],'ChecksumsPerFrame',1);
%     codeword = step(gen,Accept_Or_Reject');
%     codewordline=codeword';
% %bits=(rand(1,Nbr_Ofbits)>0.5);
% %bits=rand(1,Nbr_Ofbits)>0.5;
% %Apply convolutionnal code
%  trellis = poly2trellis(6,[77 , 45]);
%  bits = convenc(In_puts,trellis);
%  
% 
%  
% 
% %Mapping between input_01 and symbols
%  
% 
%  
% 
%        
% 
%  
% 
% vec_1Andminusone=2*bits-1; % y is a vector of 1 and -1
%  
% 
%  
% 
% Re=vec_1Andminusone(1:2:length(bits)-1); 
% Im=vec_1Andminusone(2:2:length(bits));
% Complex_sig=Re+1i*Im; % QAM Symbols 
% %1. Mapping bits to symbols
% %Transmitted_bits = Mapping_QAM_symb(Complex_sig );
% %Transmitted_bits2=Transmitted_bits';
% % 2. Generating pilots
% Adding_Pilots = adding_pilots(Nsc)';
% % 3. Adding the pilot at the beginning of the symbols
%  
% 
%  
% 
% symbols=[Adding_Pilots Complex_sig];
%  
% 
%  
% 
% % 4.Converting Serial-to-parallel by using reshape commond
% Serial_to_parallel=reshape(symbols,[Nsc,length(symbols)/Nsc]);
%  
% 
%  
% 
% % 5. Taking the IFFT of the 4
% Result_Of_ifft=ifft(Serial_to_parallel,Nsc);
%  
% 
%  
% 
% % 6.Converting Parallel to Serial by using reshape commond
% Paralle_To_Serial=reshape(Result_Of_ifft,[1,Nsc*length(symbols)/Nsc]);
%  %Remove_zeros = Removing_Paddad_Zeros(Zero_Paddad_Data)
% %7.Finding the CPs which is the last 20 values of each Nsc placed at the
% %beginning)
% Pilot = Paralle_To_Serial(1:1:128);
% Paralle_To_Serial_WithPilot= Paralle_To_Serial(129:1:end);
% for l=1:length(Paralle_To_Serial_WithPilot)/Nsc
%     CP(l,:)=Paralle_To_Serial_WithPilot(l*Nsc-19:1:l*Nsc);
% end
% %8.The CPs inserts
% Adding_Cp_WithData=Paralle_To_Serial_WithPilot;
% Adding_Cp_WithData=[CP(1,:) Adding_Cp_WithData];
% for m=1:1:length(CP(:,1))-1
%     Adding_Cp_WithData=[Adding_Cp_WithData(1:1:m*Nsc+m*20) CP(m+1,:) Adding_Cp_WithData(m*Nsc+m*20+1:1:end)];
% end 
% %9.Adding pilot again
%  Add_Pilot_Data_WithCp =[ Pilot Adding_Cp_WithData] ; 
% %10.Digital to Analog convertion 
%  
% 
%  
% 
% %Digital_To_Analog=interp1(1:1:length(Add_Pilot_Data_WithCp), Add_Pilot_Data_WithCp, 1:1/22:length(Adding_Cp_WithData));
% Digital_To_Analog=interp(Add_Pilot_Data_WithCp,20);
% %Digital_To_Analog=Digital_To_Analog(1:end-19);
% figure(1)
% plot(abs(Digital_To_Analog))
% % Modulation 
% Real_Part_Ofsig = real(Digital_To_Analog);
% Imag_Part_Ofsig = imag(Digital_To_Analog);
%  
% 
%  
% 
% %t1=0:(1/Sampling_Freq):(length(Real_Part_Ofsig)-1)/Sampling_Freq;
% t1=0:(Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq)):(length(Real_Part_Ofsig)-1)*Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq);
% % RQM=Real_Part_Ofsig.*cos(2*pi*t1*Carrier_freq);
% % RIM =-1*Imag_Part_Ofsig.*sin(2*pi*t1*Carrier_freq);
% %Transmitting_Signa= RQM+RIM;
% Carrier = exp(1i*2*pi*Carrier_freq*t1);
% Transmitting_Signa=Digital_To_Analog.*Carrier;
%  
% 
%  
% 
% % Noise = 0.01*randn(1,1000);
%  %Transmitting_Signa=[ Transmitting_Signa]; 
% Transmitting_Signa=real( Transmitting_Signa);
%  figure(2)
%  plot(Transmitting_Signa)
% title(' Trasnmitted signal');
%  
% 
% %Transmitting via microphone 
%  Sending=audioplayer(Transmitting_Signa,Sampling_Freq); 
% playblocking(Sending);
% end