function [sequence,output,energy_asf] = predict_rna_ss(seq_id,options)
% runs the fold subroutine to calculate the energy of probe folding

    if (nargin < 1)
        error('Not enough input arguments.');
    elseif (nargin < 2)
        options = '';
    end

    % Set $DATAPATH to location of thermodynamic parameter files
    setenv('DATAPATH', 'C:\Program Files\RNAstructure6.0.1\data_tables\')
    %get the sequence output necessary for this function
    [~,sequence] = read_fasta(['Probes/' seq_id '.fna']);
    % Use absolute path to desired executable
    Fold = strcat(['"C:\Program Files\RNAstructure6.0.1\exe\Fold.exe" ' ...
        , 'Probes/' seq_id '.fna', ' ', ' output/' seq_id '.ct' options]);
    
    % Run the full command
    [~,output] = system(Fold)
    
    %'fold' RNAstructure predict the lowest free energy structure and a set of sub structures

% [~,sequence] = read_fasta(['Probes/' seq_id '.fna']);
% command = ['bin/fold ' ' Probes/' seq_id '.fna' ' output/' seq_id '.ct' options];
% [~,output] = system(command);    
energy_asf = read_ct(['output/' seq_id '.ct']);

end

%for some reason only the sequence is getting outputted? 20181024