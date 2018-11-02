% function energy_asf = wrapper_energy(seq_id_list,intervals,sequence)
% overly complicated? see commented out line 24
function energy_asf = wrapper_energy(intervals,sequence)

%seq_id is the name of the probe (cell variable just as you pull it out using [~,~,seq_id] = xlsread...)
%intervals is two columns: first is starts and second is ends
%sequence is the target sequence in letters (char variable)

i = intervals(:,1);
j = intervals(:,2); 
energy_asf=zeros(length(i),1);

for m=1:length(i)
    region = [];
    for n=i(m):j(m)
    region = strcat(region,sequence(n)); 
    region = char(region);
    end
    region
probe = seqrcomplement(region); 
probe

% seq_id = char(seq_id_list(m)); overly complicated?
probe_id = num2str(m);
% fasta_builder(probe,['Probes\' seq_id '.fna'],seq_id); corrected SM
fasta_builder(probe,['Probes\' probe_id '.fna'])

[~,~,energy_asf(m)]= predict_rna_ss(probe_id,'');

end
end

