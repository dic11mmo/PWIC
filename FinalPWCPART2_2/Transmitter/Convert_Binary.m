
function Conver_To_Binary = Convert_Binary(Convert_To_10, Packet_Size) 
Data_Size= Convert_To_10(:)';
[Row,Colums] = size(Convert_To_10);
Nbr_Of_Packets= floor (length(Data_Size)/Packet_Size);
 Change_Bin_NbrPackets = de2bi(Nbr_Of_Packets,12);
 Change_Bin_lengthPacket=de2bi(Packet_Size,12);
 Change_Bin_Row =de2bi(Row,12);
 Change_Bin_Colums = de2bi(Colums,12);
  Change_Bin_Data_Size= de2bi(length(Data_Size),20);
 All_Information = [ Change_Bin_NbrPackets Change_Bin_lengthPacket  Change_Bin_Row  Change_Bin_Colums  Change_Bin_Data_Size ];
 Nbr_Zeros = Packet_Size-length(All_Information) ;
 Conver_To_Binary = [  All_Information round(randi([0,1],1, Nbr_Zeros))];