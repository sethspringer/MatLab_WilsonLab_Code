%% Return the exact p value (one tail and 2 tails), from a given z (standardized distribution)
% INPUT:
%       z value 
% OUTPUT:
%       one tail and two tails exact p value for that z
function [pvalue] = pvaluefromz_new(z)


    z = abs(z);
    F = @(x)(exp (-0.5*(x.^2))./sqrt (2*pi));
    p = integral(F, z, 100);
    
    %Commenting this out bc it is unnecessary
    %fprintf ('\nOne tail p-value : %1.6f', p);

    %fprintf ('\nTwo tail p-value : %1.6f\n' ,p*2)

    %Return two tail p-value
    pvalue = p*2; 
    
end