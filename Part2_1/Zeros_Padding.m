function Adding_Zeros = Zeros_Padding(Data_Vector , Nsc)
[Row,Colum]=size(Data_Vector); Data_Vector=Data_Vector(:).';
NbR_Of_Zeros_ToAdd_Colum =Nsc - rem(length(Data_Vector),Nsc); padd= round(rand(1,NbR_Of_Zeros_ToAdd_Colum));
Change_LEN_Data_Binary=dec2bin(length(Data_Vector),30);
Change__Row_Binnary= dec2bin(Row,30);
Change_Colum_Binary=dec2bin(Colum,30);
Vect_Of_All = [Change_LEN_Data_Binary Change__Row_Binnary Change_Colum_Binary];
Vect_Double=[];
for i=1:length(Vect_Of_All)
Vect_Double=[Vect_Double str2double(Vect_Of_All(i))]; end
NbR_Of_Zeros_ToAdd_Colum1 =Nsc-length(Vect_Of_All)-6; padd1= zeros(1,6);
Len_Dat_Row_Colums=[ Vect_Double round(rand(1, NbR_Of_Zeros_ToAdd_Colum1)) padd1];
Adding_Zeros = [Len_Dat_Row_Colums Data_Vector padd ];
end