function seqinnum = seqtonum(filename,sheet,data_range)

[~,~,sequence] = xlsread(filename,sheet,data_range);
sequence = char(sequence);
seqinnum = zeros(length(sequence),1);

for i = 1:length(sequence)
    switch sequence(i)
        case 'T'
            seqinnum(i) = 4;
        case 'U'
            seqinnum(i) = 4;
        case 'A'
            seqinnum(i) = 1;
        case 'C'
            seqinnum(i) = 2;
        case 'G'
            seqinnum(i) = 3;
        otherwise
    end
end
