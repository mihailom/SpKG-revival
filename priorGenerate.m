function [X, IN, C,basis_PT,basis_TF] = priorGenerate(T,pair,theta_0,basis_asf)
% parameters 
beta_0 = 0.39728;

% initialization
%MKM IN stands for interval
IN = [];
basis  = [];
basis_PT = [];
basis_TF = [];
N = length(T);
Cr = zeros(N);
p  = length(theta_0);

%MKM Here Yan is going through all of the possible interval combinations
%... and calculating associated energies

for k = 9:16
    for i = 1:N-k+1
        IN = [IN; i,i+k-1];
       
        [base,base_PT,base_TF] =  wrapper(i,i+k-1,T,pair);
        basis  = [basis;base];
        basis_PT = [basis_PT; base_PT];
        basis_TF = [basis_TF; base_TF];
        
    end
end
basis = basis+basis_asf;
X = repmat(basis, 1, p);

for i = 1:N
     for j = 1:N
         Cr(i,j) = exp(-beta_0*abs(i-j));
     end
end

r = 0.2;
sigma = r*theta_0;
sigma(sigma==0) = r*mean(theta_0);
C   = diag(sigma)*Cr*diag(sigma);

end