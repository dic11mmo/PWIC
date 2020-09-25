function Transmitting_Signa = TX( In_puts )
load('randomvector1.mat');
%Parameters description
Sampling_Freq = 44100;
Carrier_freq = 3000;    %carrier frequency
%Nbr_Ofbits=2^10;
Nsc = 128; %number of subcarriers in the OFDM system
Time_of_OFDMs = 58e-3;
Nbr_Of_Samples = round(Time_of_OFDMs*Sampling_Freq); % Number of samples in the OFDM information symbol time
N_Samples = round(Nbr_Of_Samples/Nsc);
 

 % to make length of data multiply of 128 
% % 
% % gen = comm.CRCGenerator([16 15 2 0],'ChecksumsPerFrame',1);
% % codeword = step(gen, In_puts');
% % codeword =codeword';
% %   In_puts = [ codeword zeros(1,880)]; 
% % %bits=(rand(1,Nbr_Ofbits)>0.5);
% % %bits=rand(1,Nbr_Ofbits)>0.5;
% %Apply convolutionnal code
In_puts = xor(randomvector1,In_puts );
 trellis = poly2trellis(6,[77 , 45]);
 bits = convenc(In_puts,trellis);
 

 

%Mapping between input_01 and symbol 
 

vec_1Andminusone=2*bits-1; % y is a vector of 1 and -1
 

 

Re=vec_1Andminusone(1:2:length(bits)-1); 
Im=vec_1Andminusone(2:2:length(bits));
Complex_sig=Re+1i*Im; % QAM Symbols 
%1. Mapping bits to symbols
%Transmitted_bits = Mapping_QAM_symb(Complex_sig );
%Transmitted_bits2=Transmitted_bits';
% 2. Generating pilots
Adding_Pilots = adding_pilots(Nsc)';
% 3. Adding the pilot at the beginning of the symbols
 

 

symbols=[Adding_Pilots Complex_sig];
 

 

% 4.Converting Serial-to-parallel by using reshape commond
Serial_to_parallel=reshape(symbols,[Nsc,length(symbols)/Nsc]);
 

 

% 5. Taking the IFFT of the 4
Result_Of_ifft=ifft(Serial_to_parallel,Nsc);
 

 

% 6.Converting Parallel to Serial by using reshape commond
Paralle_To_Serial=reshape(Result_Of_ifft,[1,Nsc*length(symbols)/Nsc]);
 %Remove_zeros = Removing_Paddad_Zeros(Zero_Paddad_Data)
%7.Finding the CPs which is the last 20 values of each Nsc placed at the
%beginning)
Pilot = Paralle_To_Serial(1:1:128);
Paralle_To_Serial_WithPilot= Paralle_To_Serial(129:1:end);
for l=1:length(Paralle_To_Serial_WithPilot)/Nsc
    CP(l,:)=Paralle_To_Serial_WithPilot(l*Nsc-19:1:l*Nsc);
end
%8.The CPs inserts
Adding_Cp_WithData=Paralle_To_Serial_WithPilot;
Adding_Cp_WithData=[CP(1,:) Adding_Cp_WithData];
for m=1:1:length(CP(:,1))-1
    Adding_Cp_WithData=[Adding_Cp_WithData(1:1:m*Nsc+m*20) CP(m+1,:) Adding_Cp_WithData(m*Nsc+m*20+1:1:end)];
end 
%9.Adding pilot again
 Add_Pilot_Data_WithCp =[ Pilot Adding_Cp_WithData] ; 
%10.Digital to Analog convertion 
 

 

%Digital_To_Analog=interp1(1:1:length(Add_Pilot_Data_WithCp), Add_Pilot_Data_WithCp, 1:1/22:length(Adding_Cp_WithData));
Digital_To_Analog=interp(Add_Pilot_Data_WithCp,20);
%Digital_To_Analog=Digital_To_Analog(1:end-19);
% figure(1)
% plot(abs(Digital_To_Analog))
% Modulation 
Real_Part_Ofsig = real(Digital_To_Analog);
Imag_Part_Ofsig = imag(Digital_To_Analog);
 

 

%t1=0:(1/Sampling_Freq):(length(Real_Part_Ofsig)-1)/Sampling_Freq;
t1=0:(Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq)):(length(Real_Part_Ofsig)-1)*Time_of_OFDMs/round(Time_of_OFDMs*Sampling_Freq);
% RQM=Real_Part_Ofsig.*cos(2*pi*t1*Carrier_freq);
% RIM =-1*Imag_Part_Ofsig.*sin(2*pi*t1*Carrier_freq);
%Transmitting_Signa= RQM+RIM;
Carrier = exp(1i*2*pi*Carrier_freq*t1);
Transmitting_Signa=Digital_To_Analog.*Carrier;
 

 

%  Noise = 0.01*randn(1,1000);
 %Transmitting_Signa=[ Transmitting_Signa]; 
Transmitting_Signa=real( Transmitting_Signa);
 figure(2)
 plot(Transmitting_Signa)
title(' Trasnmitted signal');
 

%Transmitting via microphone 
 Sending=audioplayer(Transmitting_Signa,Sampling_Freq); 
playblocking(Sending);



 end