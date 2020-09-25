% Nbr of subcarrier
function Adding_Pilots = adding_pilots(Nsc)
 

    % Obtaining the known pilots
    Adding_Pilots = zeros (Nsc,1);
    randn('state', 100);
    %rng(100,'v5normal');
    x = 2*sign(randn(1, Nsc/2));
    Adding_Pilots(1:2:end) = [x];
    

end
