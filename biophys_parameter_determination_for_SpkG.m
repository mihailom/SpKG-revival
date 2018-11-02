function [basis,X,C,theta_0] = biophys_parameter_determination_for_SpkG(excel_name,sheet_number,seq_data_range, bpp_data_range, structure_data_range,l)
% sequence data range is the number of the nucleotide and the number of its
% 
% partner -- HOW DID I GET THESE? RNA STRUCTURE? NUPACK?
% removed intervals as an input because all of the possible intervals will
% ... be generated here
% intervals =[]
% % output matrix columns: G_ast, G_tf, G_asf, thetha_bar
% for k = 9:16
%     for i = 1:l-k+1
%         intervals = [intervals; i,i+k-1];
%         
% % %         [base,base_PT,base_TF] =  wrapper(i,i+k-1,T,pair);
% % %         basis  = [basis;base];
% % %         basis_PT = [basis_PT; base_PT];
% % %         basis_TF = [basis_TF; base_TF];
%         
%     end
% end   

%test case
intervals = [1,9;35,50;66,80]
probe_specs = intervals;

    %CT_info is called "structure" in Yan's codes = [T; pair] = [num of seq
    %... nt; num of the nt that pairs with seq nt]
    CT_info = [seqtonum(excel_name,sheet_number,seq_data_range),xlsread(excel_name,sheet_number,structure_data_range)];
    
    %THETA_0 FOR PRIORS --> one bpp value for each nt in the target RNA seq
    theta_0 = xlsread(excel_name,sheet_number,bpp_data_range);
    
    %ENERGY OF PROPOSED INTERACTION
    %Yan's wrapper function did not loop through all of the intervals so I
    % used the wrapper function from the biophys model ("wrapper_MKM"), 
    % modified v slightly to call on Yan's Phi2
    % note that in Yan's codes, beta was set to 0 so I just left it like
    % that
    G_ast = wrapper_MKM(probe_specs, CT_info(:,1));
    
    %This is how we originally calculated G_ast for the biophys model,
    %using Phi (not Phi2)in the biophysical model
    %G_ast = wrapper(probe_specs,seqtonum(excel_name,sheet_number,seq_data_range));

    [~,~,sequence_as_column_vector] = xlsread(excel_name,sheet_number,seq_data_range);
    sequence_as_string = strcat(sequence_as_column_vector);
    sequence_as_string = char(sequence_as_string);
    
    %ENERGY OF TARGET REGION UNFOLDING
    % this is base_TF in Yan's code.. Gtarget region is essentially
    % the same function 20181022
    G_tf = Gtargetregion(probe_specs,CT_info,l);
    
    %ENERGY OF ASRNA UNFOLDING    
    G_asf = wrapper_energy(probe_specs,sequence_as_string');

% SUM OF BPP PER TARGET REGION
% this was for when we were summing bpp for each target region
% theta_bar =zeros(length(probe_specs(:,1)),1)
% for n=1:length(probe_specs(:,1))
% % this is theta_0 in Yan's code???  20181022  
% theta_bar(n) = sum(bpp_data(probe_specs(n,1):probe_specs(n,2)))
% end


    %probe_specs = [probe_specs, G_ast, G_tf, G_asf]
    
%start of Yan's additions
%need to figure out what theta_0 is
Cr = zeros(l);
beta_0 = 0.39728;
r = 0.2;


% Yan's "basis" variable
basis = [];
basis = G_ast + G_asf + G_tf - 1;
%WHY DOES YAN SUBTRACT 1 HERE? (from wrapper.m)
p = length(theta_0);
X = repmat(basis, 1, p);



for i = 1:l
     for j = 1:l
         Cr(i,j) = exp(-beta_0*abs(i-j));
     end
end


sigma = r*theta_0;
sigma(sigma==0) = r*mean(theta_0);
C   = diag(sigma)*Cr*diag(sigma);


end