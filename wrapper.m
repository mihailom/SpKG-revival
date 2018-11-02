function [base,base_PT,base_TF] = wrapper(i,j,T,pair)
% MKM T is a vector of nucleotides from 5' to 3' with each nucleotide 
%...listed in terms of a number (use seqtonum to make this this?)
% i,j are intervals from priorGenerate
% pair is the number of the nucleotide that nucleotide n is bound to (0 if
% ...nothing)

l = length(T);
structure = [T;pair]';
base_PT = sum(Phi2(i,j,T,0));
base_TF = Gtargetregion([i,j],structure,l);
% probe = '';
% for k = j:-1:i
%     switch T(k)
%         case 1
%             NT = 'U';
%         case 2
%             NT = 'G';
%         case 3
%             NT = 'C';
%         case 4
%             NT = 'A';
%         otherwise
%             fprintf('Invalid sequence\n');
%     end
%     probe = strcat(probe,NT);
% end
% if ~exist(['Probes/' seq_id '.fna'],'file')
%     fasta_builder(probe,['Probes/' seq_id '.fna']);
% end  
% [~,~,base_asf] = predict_rna_ss(seq_id,'');
base = base_PT+base_TF-1;
end