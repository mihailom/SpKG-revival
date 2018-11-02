function fasta_builder(probe,filename)
%builds a fasta file using a new probe sequence

EcoRI = 'GAAUUC';
downStream = 'UACCAUUCACCUCUUGGAUUUGGGUAUUAAAGAGGAGAAAGGUACCAUGAGUAAAG';
s = strcat(EcoRI,probe,downStream);
fid = fopen(filename,'wt'); 
fprintf(fid,s); 

end