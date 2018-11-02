% This function returns the amount of binding to target region [i,j] in molecule T, given 
% a probability of accesibility vector theta. We use the following encoding:
% A = 1, C = 2, G = 3, U = 4
% INPUTS:
%   i      = start of target region (between 1 and length(T), inclusive, and smaller than j)
%   j      = end of target region (between 1 and length(T), inclusive, and larger than i)
%   T      = target molecule, a vector whose entries are one of {1, 2, 3, 4},
% OUTPUTS:
%   retVal = a (1+length(T)) x 1  vector of basis function values
%   make it a row vector output % added by Yan Li 10/12/2015

function retVal = Phi2(i,j,T,beta)
    % energy parameters
    GAMMA_0 = 4.09;
    %GAMMA_0 = 0;
    GAMMA_AU = 0.45;
    GAMMA_SYM = 0.43;
    GAMMA = [ -0.93 -2.24 -2.08 -1.10; -2.11 -3.26 -2.36 -2.08; -2.35 -3.42 -3.26 -2.24; -1.33 -2.35 -2.11 -0.93];
    COMPLIMENT = [4 3 2 1];
    
    % This function extracts the sequence from an excel file
    % T = xlsread('GIintronsequence1.xlsx');
    retVal = zeros(1,1+length(T),1); % edited by Yan Li 10/12/2015

    % Base value
    retVal(1) = GAMMA_0;

    % check for self-complimentary
    isSym = 1;
    for k=0:(j-i)
        if(T(j-k) == 0)
            isSym = 0;
            break;
        end
        
        if(T(i+k) ~= COMPLIMENT(T(j-k)))
            isSym = 0;
            break;
        end
    end

    % if self-complementary, add penalty term
    if(isSym)
        retVal(1) = retVal(1) + GAMMA_SYM;
    end

    % check for AU ending penalty
    if((T(i) == 4) | (T(i) == 1))
%     if((T(i) == 4) || (T(i) == 1))
        retVal(1) = retVal(1) + GAMMA_AU;
    end
    if((T(j) == 4) | (T(j) == 1))
        retVal(1) = retVal(1) + GAMMA_AU;
    end


    % now put in energetic contribution of each base-paired stack
    for k=i:j-1
        if((T(k) == 0) | (T(k+1) == 0))
            continue;
        end
        
        retVal(k+1) = retVal(k+1) + GAMMA(T(k), T(k+1));
    end
    
    % penalize length
    %commenting this out because in Yan's "wrapper function", beta = 0
    %retVal = retVal*exp(-beta*(j-i));
    
    % finally negate to turn to maximization problem
    retVal = -retVal;
end
