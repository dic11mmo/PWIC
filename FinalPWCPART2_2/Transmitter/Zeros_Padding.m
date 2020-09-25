%%% PWC PART 2 TASK 1 
% To make the length of each massege or file be equal to multiple of Nsc We 
% are padding zeros. 
 

 

%  Data : File or packet: The matrix of packets
%   Nsc: Nbr of subcarriers
 

function Adding_Zeros = Zeros_Padding(Data_Vector , Packet_size  )
Data_Vector=Data_Vector(:).';
 NbR_Of_Zeros_ToAdd_Colum =Packet_size - rem(length(Data_Vector),Packet_size);
   padd= round(rand(1,NbR_Of_Zeros_ToAdd_Colum)); 

  Adding_Zeros = [ Data_Vector padd ];

end 
%Nsc = Packet_size;
%  [Row,Colum]=size(Data_Vector);
% Data_Vector=Data_Vector(:).';
%  
% 
%  Nbr_Of_Packets= floor (length(Data_Vector)/Packet_size);
%  NbR_Of_Zeros_ToAdd_Colum =Nsc - rem(length(Packet_size),Nsc);
%    padd= round(rand(1,NbR_Of_Zeros_ToAdd_Colum)); 
%    
% 
%    length_of_Packet_Binary =dec2bin(Packet_size ,30); 
%     Nbr_Of_Packets_Binnary= dec2bin(Nbr_Of_Packets,30);
%       %LEN_Data_Binary=dec2bin(length(Data_Vector),30);
%         Row_Binnary= dec2bin(Row,30);
%            Colum_Binary=dec2bin(Colum,30);
%     Vect_Of_All = [length_of_Packet_Binary Nbr_Of_Packets_Binnary Row_Binnary Colum_Binary];
%  
% 
%     Vect_Double=[];
% for i=1:length(Vect_Of_All)
%     Vect_Double=[Vect_Double str2double(Vect_Of_All(i))];
% end 
%  
% 
% NbR_Of_Zeros_ToAdd_Colum1 =Nsc-length(Vect_Of_All)-6;              
%                  padd1= zeros(1,6);
%                  Len_Dat_Row_Colums=[ Vect_Double round(rand(1, NbR_Of_Zeros_ToAdd_Colum1)) padd1];
%                  Adding_Zeros = [Len_Dat_Row_Colums Data_Vector padd ];

