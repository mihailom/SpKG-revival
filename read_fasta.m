function [description,sequence] = read_fasta(filename)
%grabs sequence from fasta file

fid = fopen(filename,'r');
description = fgetl(fid);
tline = fgetl(fid);
sequence = '';
while ischar(tline)
    sequence = [sequence tline];
    tline = fgetl(fid);
end

fclose(fid);
end