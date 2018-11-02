function [Gibbs] = wrapper_MKM(intervals, T) 
%intervals are provided as a double variable with two columns column 1 is starts and column 2 is ends of regions
%T is the sequence of the target RNA as a column (only one column in numbers -you can use seqtonum.m function to do this)
Gibbs = zeros(length(intervals),1);
    
    for n = 1 : length(intervals)
        i = intervals (n,1);%- 21; %take the 21 out when using other molecules
        j = intervals (n,2);%- 21; %take the 21 out when using other molecules 
        retVal = Phi2(i,j,T,0);
        % Gibbs Free Energy of binding to target area
        Gibbs(n,1) = sum (retVal);
       
    end
    
end
    