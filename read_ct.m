function [energy_asf] = read_ct(filename)
%grabs the energy value from a ct file

fileID = fopen(filename);
C = textscan(fileID, '%s');
%celldisp(C);
energy_asf = str2double(C{1,1}(4));
fclose(fileID);
end

