
clear all
close all
clc

Time2 = 0;
 tic;
%     randomvector1 = randi([0,1],1,1024);
%  save randomvector

Recording_Time=1.7;

Sampling_Freq = 44100;
Packet_Size = 2^10;
Packet_Size = Packet_Size-31;

image = imread('forMatlab.jpg');
Convert_To_10= im2bw(image); 
Data_Vector= Convert_To_10(:)';
  Nbr_Of_Packets= ceil(length(Data_Vector)/Packet_Size);%rounding number of packets to positive infinity as we dont want to lose data%

 Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
 Padd_Zeros_With_Data = Zeros_Padding(Convert_To_10, Packet_Size );
 
 Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size);
  Matrix_Of_Packets(1,:) =  Conver_To_Binary ;

 a=2;
 % Starts inserting data in the second row of matrix
for l=1:993:length(Padd_Zeros_With_Data)-992
     Matrix_Of_Packets(a,:)= Padd_Zeros_With_Data( l*1:l+(Packet_Size-1));
     
     a=a+1;
end

 while Time2 < 1
    Time2 = toc
end

 k = 1; 
 i=1;
 while i<=Nbr_Of_Packets
  Time = 0;
    tic; 
    if i==1
         disp(' The first packet is sending')
         disp(i)
      Packet= Matrix_Of_Packets(i,:); 
 Header= de2bi(k,15);
  Len_Of_Packets = [Header Packet]; % The length of
     Data_To_transmit=TX(Len_Of_Packets);
     i=i+1;
%      pause(0.5)
     end
     
disp('record start')
recorder=audiorecorder(Sampling_Freq,16,1);
 recordblocking(recorder,Recording_Time)
R=getaudiodata(recorder); 
disp('record stop')
 
     Received_Data =  RX(R);
  Ack_Or_Nack= mean (Received_Data);
   if  Ack_Or_Nack >0.4 % Assum it Nack if it is > 0.4  
       i=i-1;
 disp(' Last packet is sending')
         disp(i)
%          k=k-1;
   Packet= Matrix_Of_Packets(i,:);
 Header= de2bi(k,15);
  Len_Of_Packets = [Header Packet]; % The length of
        pause(0.5)
    Data_To_transmit=TX(Len_Of_Packets); 
     disp('The packet has been send ')
    i=i+1;
   elseif Ack_Or_Nack <=0.4 
   k = k+1 ;
      i= i+1 ;
  Packet= Matrix_Of_Packets(i,:);
 Header= de2bi(k,15);
  Len_Of_Packets = [Header Packet]; 
     Data_To_transmit=TX(Len_Of_Packets); 
     disp('Send packet nbr ')
     disp(i)
   end
 while Time < 2.4
        Time = toc;
    end

 end 






% clear all
% close all
% clc
% 
% Time2 = 0;
% tic;
% % load('randomvector.mat')
% Recording_Time=2;
% Sampling_Freq = 44100;
% Packet_Size = 2^10;
% Packet_Size = Packet_Size-31;
% image = imread('forMatlab.jpg');
% Convert_To_10= im2bw(image); 
% Data_Vector= Convert_To_10(:).';
%   Nbr_Of_Packets= floor(length(Data_Vector)/Packet_Size);%rounding number of packets to positive infinity as we dont want to lose data%
% 
%   Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  Padd_Zeros_With_Data = Zeros_Padding(Convert_To_10, Packet_Size );
%  
%   Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size);
%   % Adding information about nbr rows, Colums length of paket .....
%   % In the  first row of Matrix_Of_Packets as follows 
%   Matrix_Of_Packets(1,:) =  Conver_To_Binary ;
%  a=2;
%  % Starts inserting data in the second row of matrix
% for l=1:993:length(Padd_Zeros_With_Data)-992
%      Matrix_Of_Packets(a,:)= Padd_Zeros_With_Data( l*1:l+(Packet_Size-1));
%      
%      a=a+1;
% end
% 
% % Matrix_Of_Packets2= transp(Matrix_Of_Packets);
% % Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
% % Rem_Zer= length(Matrix_Of_Packets22)-length(Data_Vector);
% % Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-Rem_Zer);
% % Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[262,193]);
% % imshow( Matrix_Of_Packets4)
% 
% 
% 
%  while Time2 < 1
%     Time2 = toc
%  end
%  
% 
%  i=1;
% %  Nbr_Of_Packets
%  while i<=1
%   Time = 0;
%     tic; 
%     if i==1
%          disp(' First packet is sending')
%          disp(i)
%          Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(i,15);
%   Len_Of_Packets = [Header Packet]; % The length of
%      Data_To_transmit=TX(Len_Of_Packets);
%      i=i+1;
% 
%      end
%      
% disp('record start')
% recorder=audiorecorder(Sampling_Freq,16,1);
%  recordblocking(recorder,Recording_Time)
% R=getaudiodata(recorder); 
% disp('record stop')
%  
%      Received_Data =  RX(R);
%   Ack_Or_Nack= mean (Received_Data);
%    if  Ack_Or_Nack >0.4 % Assum it Nack   
%        i=i-1;
%  disp(' Received Nack last packet is sending')
%          disp(i)
% 
%    Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(i,15);
%   Len_Of_Packets = [Header Packet]; % The length of
%       pause(0.5)
%     Data_To_transmit=TX(Len_Of_Packets); 
%   
%      disp('The packet has been send ')
%     i=i+1;
%    elseif Ack_Or_Nack <=0.2
%    i = i+1 ;
%       i= i+1 ;
%   Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(i,15);
%   Len_Of_Packets = [Header Packet]; 
%      Data_To_transmit=TX(Len_Of_Packets); 
%      disp('Send packet nbr ')
%      disp(i)
%    end
%  while Time < 2.4
%         Time = toc
%     end
% 
%  end 


















































% clear all
% close all
% clc
% 
% Time2 = 0;
% tic;
% 
% % load('randomvector.mat')
% 
% Recording_Time=1.6;
% 
% Sampling_Freq = 44100;
% % Nsc = 128;
% Packet_Size = 2^10;
% Packet_Size = Packet_Size-31;
% %  image = imread('images.jpg');
% image = imread('forMatlab.jpg');
% Convert_To_10= im2bw(image); 
% Data_Vector= Convert_To_10(:).';
%   Nbr_Of_Packets= floor(length(Data_Vector)/Packet_Size);%rounding number of packets to positive infinity as we dont want to lose data%
% %    Nbr_Of_Packets=5;
%   Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  Padd_Zeros_With_Data = Zeros_Padding(Convert_To_10, Packet_Size );
% 
%   %Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  %.1
%   Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size);
% %.1 end
% % Conver_To_Binary = Convert_Binary_Witdec2bin(Convert_To_10, Packet_Size);
% %   Matrix_Of_Packets(1,:) =  Conver_To_Binary ;
%  a=1;
% for l=1:993:length(Padd_Zeros_With_Data)-992
%      Matrix_Of_Packets(a,:)= Padd_Zeros_With_Data( l*1:l+(Packet_Size-1));
%      
%      a=a+1;
% end
% 
% %%% If i convert Matrix_Of_Packets to a matrix of 262x193 like Matrix_Of_Packets1 = Matrix_Of_Packets (1:end-77) 
% % here 77 is the zeros which is added by zero_padding function and using the
% % imshow( matrisof262x193) i didnt get the picture back  solved as followed
% 
% % Matrix_Of_Packets2= transp(Matrix_Of_Packets);
% % Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
% % Rem_Zer= length(Matrix_Of_Packets22)-length(Data_Vector);
% % Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-Rem_Zer);
% % Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[262,193]);
% % imshow( Matrix_Of_Packets4)
% 
% 
% 
%  while Time2 < 1
%     Time2 = toc
% end
% 
%  k = 1; 
%  i=1;
%  while i<=Nbr_Of_Packets
%   Time = 0;
%     tic; 
% %     if i==1
%          disp(' First packet is sending')
%          disp(i)
%          Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(k,15);
%   Len_Of_Packets = [Header Packet]; % The length of
%      Data_To_transmit=TX(Len_Of_Packets);
%      i=i+1;
%       k=k+1;
%       pause(1)
% %      end
% %      
% % % disp('record start')
% % % recorder=audiorecorder(Sampling_Freq,16,1);
% % %  recordblocking(recorder,Recording_Time)
% % % R=getaudiodata(recorder); 
% % % disp('record stop')
% %  
% %      Received_Data =  RX(R);
% %   Ack_Or_Nack= mean (Received_Data);
% %    if  Ack_Or_Nack >0.4 % Assum it Nack   
% %        i=i-1;
% % %          k=k-1;
% %  disp(' Received Nack last packet is sending')
% %          disp(i)
% % %          k=k-1;
% %    Packet= Matrix_Of_Packets(i,:);
% %  Header= de2bi(k,15);
% %   Len_Of_Packets = [Header Packet]; % The length of
% % %         pause(0.2)
% %     Data_To_transmit=TX(Len_Of_Packets); 
% %      disp('The packet has been send ')
% %     i=i+1;
% %    elseif Ack_Or_Nack <=0.4 
% %    k = k+1 ;
% % %       i= i+1 ;
% %   Packet= Matrix_Of_Packets(i,:);
% %  Header= de2bi(k,15);
% %   Len_Of_Packets = [Header Packet]; 
% %      Data_To_transmit=TX(Len_Of_Packets); 
% %      disp('Send packet nbr ')
% %      disp(i)
% %    end
% %  while Time < 2.4
% %         Time = toc
% %     end
% 
%  end 





































% 
% clear all
% close all
% clc
% 
% Time2 = 0;
% tic;
% 
% % load('randomvector.mat')
% 
% Recording_Time=1.6;
% 
% Sampling_Freq = 44100;
% % Nsc = 128;
% Packet_Size = 2^10;
% Packet_Size = Packet_Size-31;
% %  image = imread('images.jpg');
% image = imread('forMatlab.jpg');
% Convert_To_10= im2bw(image); 
% Data_Vector= Convert_To_10(:).';
%   Nbr_Of_Packets= floor(length(Data_Vector)/Packet_Size);%rounding number of packets to positive infinity as we dont want to lose data%
% %    Nbr_Of_Packets=5;
%   Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  Padd_Zeros_With_Data = Zeros_Padding(Convert_To_10, Packet_Size );
% 
%   %Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  %.1
%   Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size);
% %.1 end
% % Conver_To_Binary = Convert_Binary_Witdec2bin(Convert_To_10, Packet_Size);
%   Matrix_Of_Packets(1,:) =  Conver_To_Binary ;
%  a=2;
% for l=1:993:length(Padd_Zeros_With_Data)-992
%      Matrix_Of_Packets(a,:)= Padd_Zeros_With_Data( l*1:l+(Packet_Size-1));
%      
%      a=a+1;
% end
% 
% %%% If i convert Matrix_Of_Packets to a matrix of 262x193 like Matrix_Of_Packets1 = Matrix_Of_Packets (1:end-77) 
% % here 77 is the zeros which is added by zero_padding function and using the
% % imshow( matrisof262x193) i didnt get the picture back  solved as followed
% 
% % Matrix_Of_Packets2= transp(Matrix_Of_Packets);
% % Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
% % Rem_Zer= length(Matrix_Of_Packets22)-length(Data_Vector);
% % Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-Rem_Zer);
% % Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[262,193]);
% % imshow( Matrix_Of_Packets4)
% 
% 
% 
%  while Time2 < 1
%     Time2 = toc
% end
% 
%  k = 1; 
%  i=1;
%  while i<=Nbr_Of_Packets
%   Time = 0;
%     tic; 
%     if i==1
%          disp(' First packet is sending')
%          disp(i)
%          Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(k,15);
%   Len_Of_Packets = [Header Packet]; % The length of
%      Data_To_transmit=TX(Len_Of_Packets);
%      i=i+1;
% %      k=k+1;
% %      pause(0.5)
%      end
%      
% disp('record start')
% recorder=audiorecorder(Sampling_Freq,16,1);
%  recordblocking(recorder,Recording_Time)
% R=getaudiodata(recorder); 
% disp('record stop')
%  
%      Received_Data =  RX(R);
%   Ack_Or_Nack= mean (Received_Data);
%    if  Ack_Or_Nack >0.4 % Assum it Nack   
%        i=i-1;
% %          k=k-1;
%  disp(' Received Nack last packet is sending')
%          disp(i)
% %          k=k-1;
%    Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(k,15);
%   Len_Of_Packets = [Header Packet]; % The length of
% %         pause(0.2)
%     Data_To_transmit=TX(Len_Of_Packets); 
%      disp('The packet has been send ')
%     i=i+1;
%    elseif Ack_Or_Nack <=0.4 
%    k = k+1 ;
%       i= i+1 ;
%   Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(k,15);
%   Len_Of_Packets = [Header Packet]; 
%      Data_To_transmit=TX(Len_Of_Packets); 
%      disp('Send packet nbr ')
%      disp(i)
%    end
%  while Time < 2.4
%         Time = toc
%     end
% 
%  end 























% clear all
% close all
% clc
% 
% Time2 = 0;
% tic;
% 
% % load('randomvector.mat')
% 
% Recording_Time=1.6;
% 
% Sampling_Freq = 44100;
% % Nsc = 128;
% Packet_Size = 2^10;
% Packet_Size = Packet_Size-31;
% %  image = imread('images.jpg');
% image = imread('forMatlab.jpg');
% Convert_To_10= im2bw(image); 
% Data_Vector= Convert_To_10(:).';
%   Nbr_Of_Packets= floor(length(Data_Vector)/Packet_Size);%rounding number of packets to positive infinity as we dont want to lose data%
% %    Nbr_Of_Packets=5;
%   Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  Padd_Zeros_With_Data = Zeros_Padding(Convert_To_10, Packet_Size );
% 
%   %Matrix_Of_Packets = zeros(Nbr_Of_Packets,Packet_Size);
%  %.1
% %  Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size ,Nbr_Of_Packets);
% %.1 end
% % Conver_To_Binary = Convert_Binary_Witdec2bin(Convert_To_10, Packet_Size);
% %   Matrix_Of_Packets(1,:) =  Conver_To_Binary ;
%  a=1;
% for l=1:993:length(Padd_Zeros_With_Data)-992
%      Matrix_Of_Packets(a,:)= Padd_Zeros_With_Data( l*1:l+(Packet_Size-1));
%      
%      a=a+1;
% end
% 
% %%% If i convert Matrix_Of_Packets to a matrix of 262x193 like Matrix_Of_Packets1 = Matrix_Of_Packets (1:end-77) 
% % here 77 is the zeros which is added by zero_padding function and using the
% % imshow( matrisof262x193) i didnt get the picture back  solved as followed
% 
% % Matrix_Of_Packets2= transp(Matrix_Of_Packets);
% % Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
% % Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-77);
% % Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[262,193]);
% % imshow( Matrix_Of_Packets4)
% 
% 
% 
%  while Time2 < 0.9
%     Time2 = toc
% end
% 
%  k = 1; 
%  i=1;
%  while i<=Nbr_Of_Packets
%   Time = 0;
%     tic; 
%     if i==1
%          disp(' First packet is sending')
%          disp(i)
%          Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size ,Nbr_Of_Packets);
%  Header= de2bi(k,15);
%   Len_Of_Packets = [Header Conver_To_Binary]; % The length of
%      Data_To_transmit=TX(Len_Of_Packets);
%      i=i+1;
% %      pause(0.5)
%      end
%      
% disp('record start')
% recorder=audiorecorder(Sampling_Freq,16,1);
%  recordblocking(recorder,Recording_Time)
% R=getaudiodata(recorder); 
% disp('record stop')
%  
%      Received_Data =  RX(R);
%   Ack_Or_Nack= mean (Received_Data);
%    if  Ack_Or_Nack >0.2 % Assum it Nack if it is > 0.3  
%        i=i-1;
%  disp(' Received Nack last packet is sending')
%          disp(i)
% %          k=k-1;
%    Packet= Matrix_Of_Packets(i,:);
%  Header= de2bi(k,15);
%   Len_Of_Packets = [Header Packet]; % The length of
%         pause(0.5)
%     Data_To_transmit=TX(Len_Of_Packets); 
%      disp('The packet has been send ')
%     i=i+1;
%    elseif Ack_Or_Nack <=0.2 
%    k = k+1 ;
%       i= i+1 ;
% %   Packet= Matrix_Of_Packets(i,:);
% %  Header= de2bi(k,15);
% %   Len_Of_Packets = [Header Packet]; 
%      Data_To_transmit=TX(Len_Of_Packets); 
%      disp('Send packet nbr ')
%      disp(i)
%    end
%  while Time < 2.4
%         Time = toc
%     end
% 
%  end


 
 





















 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
% %  Data_Recevd=zeros(Nbr_Of_Packets ,Packet_Size);
% %  Ack_Or_Nack=zeros(Nbr_Of_Packets,32);
%  Change_Bin_NbrPackets = de2bi(Nbr_Of_Packets,30);
%  Change_Bin_lengthPacket=de2bi(Packet_Size,30);
%  Change_Bin_Row =de2bi(Row,30);
%  Change_Bin_Colums = de2bi(Colums,30);
%  All_Information = [ Change_Bin_NbrPackets Change_Bin_lengthPacket  Change_Bin_Row  Change_Bin_Colums ];
%  Nbr_Zeros = Packet_Size-length(All_Information) ;
%    All_Information = [  All_Information round(randi([0,1],1, Nbr_Zeros))];













