

clear;
close all;
clc;
Time2 = 0;
tic;
%  load('randomvector1.mat');
Recording_Time = 2;
recorder=audiorecorder(44100,16,1);
 Length_Of_The_Packet= 10; % Real information bit 
  number_of_packets = 10 ;
  ROWS= 10;
  Colums= 10;
Last_Packets_Header =0;

Nbr_Nack = 0;
%    Received_Data = zeros( number_of_packets, Length_Of_The_Packet);
%%%% Waiting one min
while Time2 < 1
    
    Time2 = toc
    
end
 i=1;   
while i<=number_of_packets
  Time = 0;
    tic; 
   disp('Start recording')
   recordblocking(recorder,Recording_Time); 
   R=getaudiodata(recorder);
 disp('recording stop')
 [CRC,Received_packet, Accept_Or_Reject ] = RX(R );
  Header = Received_packet(1:15);
  Header = bi2de(Header);
    if  CRC==0 && Header==1
       Received_packet = Received_packet(16:end);

        Nbr_packet =  Received_packet(1:12);
        Nbr_Packets = bi2de(Nbr_packet);
         number_of_packets = Nbr_Packets;
          size_Paket =  Received_packet(13:24);
           Size_Packet = bi2de(size_Paket);
           Length_Of_The_Packet = Size_Packet;
           Nbr_row =  Received_packet(25:36);
           Nbr_Rows1 = bi2de( Nbr_row);
           ROWS= Nbr_Rows1;
            Nbr_colums =  Received_packet(37:48);
           Nbr_Colums = bi2de(Nbr_colums);
            Colums=  Nbr_Colums;
           Data_Size = Received_packet(49:68);
            Data_Size= bi2de( Data_Size);

           Last_Packets_Header = Last_Packets_Header+1; 
           i=1+1;
          Data_To_transmit=TX(Accept_Or_Reject);
              disp(' Ack send the next packet ');
        Received_Data = zeros( number_of_packets, Length_Of_The_Packet);      
   elseif   CRC==0 && Header >Last_Packets_Header  
        Received_packet= Received_packet(16:end);
    Received_Data(i,:)=Received_packet;% Save received data in data 

     Last_Packets_Header = Last_Packets_Header+1;
    i=i+1;
            
     Data_To_transmit=TX(Accept_Or_Reject);
    disp('  CRC is zero Ack send the next packet ');
   else
        Data_To_transmit=TX(Accept_Or_Reject);
        disp(' Nack send the packet again');
        Nbr_Nack = Nbr_Nack+1;
       pause(0.1)
   end
 while Time < 2.4
        Time = toc;
    end
    while Nbr_Nack>=7 
        disp('Sorry check you connection  ') 
    
    end
end
 Matrix_Of_Packets2= transp(Received_Data);
Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
 Remove_Padded_Zeros = length(Matrix_Of_Packets22) -length(Data_Size);
Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-Remove_Padded_Zeros);
Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[ROWS, Nbr_Colums]);
imshow( Matrix_Of_Packets4)






























% clear;
% close all;
% clc;
% Time2 = 0;
% tic;
% %  load('randomvector1.mat');
% % load('Received_Data');
% % save Received_Data
% Recording_Time = 2;
% recorder=audiorecorder(44100,16,1);
%  Length_Of_The_Packet= 993; % Real information bit 
%   number_of_packets = 50 ;
% Nbr_Rows=262;
% Nbr_Colums=193;
% Data_Size = 50566;
% % 993 zeros + 15 Header sending to TX and Tx adds 16 CRC = 1024
% 
% %  Received_Data = zeros( number_of_packets, Length_Of_The_Packet);
% %  save Received_Data 
% Nbr_Nack = 0;
% %    Received_Data = zeros( number_of_packets, Length_Of_The_Packet);
% %%%% Waiting one min
% while Time2 < 1
%     
%     Time2 = toc
%     
% end
% Last_Packets_Header =16;
%  i=1;   
% while i<=number_of_packets
%   Time = 0;
%     tic; 
%    disp('Start recording')
%    recordblocking(recorder,Recording_Time); 
%    R=getaudiodata(recorder);
%  disp('recording stop')
%  [CRC,Received_packet, Accept_Or_Reject ] = RX(R );
%   Header = Received_packet(1:15);
%   Header = bi2de(Header);
%   disp(Header)
% %     if   Header==1
% % %         CRC==0 && Header==1
% %        Received_packet = Received_packet(16:end);
% % % [number_of_packets Length_Of_The_Packet Nbr_Rows Nbr_Colums  Data_Size]= Convert_to_decimal(Received_packet)
% %        % 1 om for Convert_to_binry function 
% %         Nbr_packet =  Received_packet(1:12);
% %         Nbr_Packets = bi2de(Nbr_packet);
% %          number_of_packets = Nbr_Packets;
% %           size_Paket =  Received_packet(13:25);
% %            Size_Packet = bi2de(size_Paket);
% %            Length_Of_The_Packet = Size_Packet;
% %            Nbr_row =  Received_packet(26:38);
% %            Nbr_Rows1 = bi2de( Nbr_row);
% %             Nbr_colums =  Received_packet(39:51);
% %            Nbr_Colums = bi2de(Nbr_colums);
% %            Data_Size = Received_packet(52:72);
% %             Data_Size= bi2de( Data_Size);
% % 
% %            Last_Packets_Header = Last_Packets_Header+1; 
% %            i=1+1;
% % %           Data_To_transmit=TX(Accept_Or_Reject);
% %               disp(' Ack send the next packet ');
% % %         Received_Data = zeros( number_of_packets, Length_Of_The_Packet);      
%    if   Header ==Last_Packets_Header+1
% %          CRC==0 && Header >Last_Packets_Header
%         Received_packet= Received_packet(16:end);
%   
%     Received_Data(i,:)=Received_packet;% Save received data in data 
%    save save Received_Data 
%      Last_Packets_Header = Last_Packets_Header+1;
%     i=i+1;
%             
% %      Data_To_transmit=TX(Accept_Or_Reject);
%     disp('  CRC is zero Ack send the next packet ');
% %    else
% %         Data_To_transmit=TX(Accept_Or_Reject);
% %         disp(' Nack send the packet again');
% %         Nbr_Nack = Nbr_Nack+1;
% %        pause(0.3)
%    end
% %  while Time < 2.4
% %         Time = toc;
% %     end
%     while Nbr_Nack>=4 
%         disp('Sorry check you connection  ') 
%     
%     end
% end
% 
%  Matrix_Of_Packets2= transp(Received_Data);
% Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
%  Remove_Padded_Zeros = length(Matrix_Of_Packets22) -length(Data_Size);
% Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-Remove_Padded_Zeros);
% Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[262,193]);
% imshow( Matrix_Of_Packets4)
% 
% % Remove_Padded_Zeros = length(Received_Data) -length(Data_Size);
% % Received_Data = Received_Data(1:end-Remove_Padded_Zeros);              
% % Convert_Matrix= reshape(Received_Data,[Nbr_Rows,Nbr_Colums]);   
% %  figure(1)
% %   imshow(Received_Data)  
















% 
% clear;
% close all;
% clc;
% Time2 = 0;
% tic;
% %  load('randomvector1.mat');
% Recording_Time = 2;
% recorder=audiorecorder(44100,16,1);
%  Length_Of_The_Packet= 9; % Real information bit 
%   number_of_packets = 1 ;
% % Nbr_Rows=262;
% % Nbr_Colums=193;
% % Data_Size = 50566;
% % 993 zeros + 15 Header sending to TX and Tx adds 16 CRC = 1024
% Last_Packets_Header =0;
% 
% Nbr_Nack = 0;
% %    Received_Data = zeros( number_of_packets, Length_Of_The_Packet);
% %%%% Waiting one min
% while Time2 < 1
%     
%     Time2 = toc
%     
% end
%  i=1;   
% while i<=number_of_packets
%   Time = 0;
%     tic; 
%    disp('Start recording')
%    recordblocking(recorder,Recording_Time); 
%    R=getaudiodata(recorder);
%  disp('recording stop')
%  [CRC,Received_packet, Accept_Or_Reject ] = RX(R );
%   Header = Received_packet(1:15);
%   Header = bi2de(Header);
%     if   CRC==0 && Header==1
%        Received_packet = Received_packet(16:end);
% % [number_of_packets Length_Of_The_Packet Nbr_Rows Nbr_Colums  Data_Size]= Convert_to_decimal(Received_packet)
%        % 1 om for Convert_to_binry function 
%         Nbr_packet =  Received_packet(1:12);
%         Nbr_Packets = bi2de(Nbr_packet);
%          number_of_packets = Nbr_Packets;
%           size_Paket =  Received_packet(13:25);
%            Size_Packet = bi2de(size_Paket);
%            Length_Of_The_Packet = Size_Packet;
%            Nbr_row =  Received_packet(26:38);
%            Nbr_Rows1 = bi2de( Nbr_row);
%             Nbr_colums =  Received_packet(39:51);
%            Nbr_Colums = bi2de(Nbr_colums);
%            Data_Size = Received_packet(52:72);
%             Data_Size= bi2de( Data_Size);
% 
%            Last_Packets_Header = Last_Packets_Header+1; 
%            i=1+1;
%           Data_To_transmit=TX(Accept_Or_Reject);
%               disp(' Ack send the next packet ');
%         Received_Data = zeros( number_of_packets, Length_Of_The_Packet);      
%    elseif  CRC==0 && Header >Last_Packets_Header
%          
%         Received_packet= Received_packet(16:end);
%     Received_Data(i,:)=Received_packet;% Save received data in data 
% 
%      Last_Packets_Header = Last_Packets_Header+1;
%     i=i+1;
%             
%      Data_To_transmit=TX(Accept_Or_Reject);
%     disp('  CRC is zero Ack send the next packet ');
%    else
%         Data_To_transmit=TX(Accept_Or_Reject);
%         disp(' Nack send the packet again');
%         Nbr_Nack = Nbr_Nack+1;
%        pause(0.3)
%    end
%  while Time < 2.4
%         Time = toc;
%     end
%     while Nbr_Nack>=4 
%         disp('Sorry check you connection  ') 
%     
%     end
% end
% 
%  Matrix_Of_Packets2= transp(Received_Data);
% Matrix_Of_Packets22= Matrix_Of_Packets2(:)';
%  Remove_Padded_Zeros = length(Matrix_Of_Packets22) -length(Data_Size);
% Matrix_Of_Packets3 = Matrix_Of_Packets22(1:end-Remove_Padded_Zeros);
% Matrix_Of_Packets4 = reshape(Matrix_Of_Packets3 ,[262,193]);
% imshow( Matrix_Of_Packets4)
% 
% % Remove_Padded_Zeros = length(Received_Data) -length(Data_Size);
% % Received_Data = Received_Data(1:end-Remove_Padded_Zeros);              
% % Convert_Matrix= reshape(Received_Data,[Nbr_Rows,Nbr_Colums]);   
% %  figure(1)
% %   imshow(Received_Data) 















