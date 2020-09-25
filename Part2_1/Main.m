clear all close all clc
tic;
%Parameters description 
Nsc=128;
% %1.
Messag_To_Send= ('For all stochastic variables: The mean of a sum equals the sum of the means. of a sum equals the sum of the mean......1 For all stochastic variables: The mean of a sum equals the sum of the means. of a sum equals the sum of the mean......2 For all stochastic variables: The mean of a sum equals the sum of the means. of a sum equals the sum of the mean......3 For all stochastic variables: The mean of a sum equals the sum of the means. of a sum equals the sum of the mean......4 ');
%
Change_To_bin=dec2bin(Messag_To_Send);
Change_To_bin=Change_To_bin.';
Change_To_bin=Change_To_bin(:);
%
input_char=[];
for i=1:length(Change_To_bin)
input_char=[input_char Change_To_bin(i)];
end
Convert_To_01=num2str(input_char)-48; %
% trellis = poly2trellis(6,[77 , 45]);
% bits = convenc(Convert_To_01,trellis); %
Padd_Zeros_With_Data = Zeros_Padding( Convert_To_01 , Nsc); %1.
% vetc2=[];
% for k =0 :(length(Padd_Zeros_With_Data )/128)-2
%
%
%
% end %
%vect= Padd_Zeros_With_Data(k*Nsc+1:k*Nsc+2*Nsc); %vetc2=[vetc2 vect ];
% Data_To_transmit=TX(vect);
Data_To_transmit=TX(Padd_Zeros_With_Data);
%% fOR iMAGE %%%% 
% image = imread('B_Whit2.jpg');
% Convert_To_10= im2bw(image); %Reducing the size of the picture by converting to gray scale %Black_Whit_Imag = dither(image); % Convert to black and white to convert the bixels values to binary values
% Data_Vector=Convert_To_10(:)';
% Padd_Zeros_With_Data = Zeros_Padding(Convert_To_10 , Nsc); Data_To_transmit=TX(Padd_Zeros_With_Data );
% 

