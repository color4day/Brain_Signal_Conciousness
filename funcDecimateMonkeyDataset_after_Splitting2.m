function funcDecimateMonkeyDataset_after_Splitting2(filenameIN, type_trial, Decim, PAR)
%  funcDecimateMonkeyDataset_after_splitting2( filenameIN, type_trial, Decim, PAR)
% 
% Example of INPUTS
% - filenameIN   = 'M_20180530-shape-rectangle-generated-on-05-04-2020 16-32.mat';
% - type_trial   = 'shape-rectangle'; % which MUST be contained in finenameIN;
% - Decim        = Decimation Factor;
%
% - [param]:    Optional set of parameters to activate filtering 
%               [Default no filtering] (see the code !!)
%
%  Script that DECIMATES and possibly FILTERs the Monkey Datasets (using'par')
%  AFTER it has been split in 6 different Datasets (each one for a specific TRIAL-TYPE)
%
% It makes use of the functions
%   - cell2matrix.m
%   - matrix2cellarray.m
% 
% to go back and forth from data stored into cells to data stored in 
% matrices. In performs decimation and/or filtering by
%
%  X = decimate_and_filtering_Matlab( Y, Fs, Decim, par );
%  
% to decimate and/or filter Signals Y by Matlab functions
%
%  - Notch-filter out the AC at 60 Hz 
%  - Low-pass filtering the signals
%  - Band-pass filtering the signals
%  All the filters are designd by Matlab "filterdesign" and 
%  implemented as zero-phase Filters by "filtfilt"
%  
%  Possibility to activate plots of signals and spectrum
%  See help of decimate_and_filtering_Matlab.m
%
% It creates an output file with the same identical structure and info
% of the original input files. The only difference is that the signals are
% (filtered and) decimated. The time-stamps are those of the
% original file. 
% Also the indexes are the same of the original sampling frequency
%

% Fs = 40000;     % Sampling Frequency
% Decim = 100;    % Decimation Factor ( 1 = No decimation)

window_name ={'baseline_sig'; 'abs_sig'; 'delay_1_sig' ; 'conc_sig' ; 'delay_2_sig' ;'choice_sig'; 'trail_sig'};
load(filenameIN)
if exist('PAR','var')
  par = PAR; % To override par if already contained in filenameIN
end
clear PAR;
Fs = fs;

if ((nargin < 4) && (~exist('param','var')))
    par.AC_removal = 0; % 1 to Remove AC by Notch-filtering at 60 Hz
    
    par.LP_Filter  = 0; % 1 to Activate Low-Pass filter with MANDATORY
    par.LP_Freqs   = [40,42]; % [Fpass, Fstop] frequencies for the filter
    
    par.BP_Filter  = 0; % 1 to Activate Band-Pass filter with MANDATORY
    par.BP_Freqs   = [28,30,70,72]; % [Fstop1 , Fpass1 , Fpass2, Fstop2]
    
    par.HP_Filter = 0;         % 1 to Activate High-Pass filter with MANDATORY
    par.HP_Freqs  = [0.4,0.5]; % [Fstop,Fpass] frequencies for the filter
    % Elliptic Filter designed specifically to High-Pass
    % 0.5 Hz, Fs = 256Hz, 40dB stopband with "filtfilt"
    % See details in the function
    
    par.PlotSignalsON = 0;     % 1 to plot original and decimated signals for comparison
    par.PlotSignalsIndex =[1,21,33,47,72,86,97,111]; % Index numbers of signals and spectrum to be plotted
    % on a 'subplot' figure
    % (Only the first 10 will be actually plotted);
    
    par.MaxSamplesPlot = inf ;            % Max # of original time-samples to be plotted (inf = ALL)
    par.MaxFreqPlot = Fs/(2*Decim); % Max frequency shown in the plot
end


for win = 1:length(window_name)
    Y       =  eval( window_name{win}); % Cell Array {signal,Trial}
    par.saveON = 0;          % 1 to save the decimated/filtered matrix to a file
    [Nch,Ntrial] = size(Y);
    
    Xcell = cell(Nch,Ntrial); % initialize to empty all the rows of the Decimated Cell-Array
    
    for trial=1:Ntrial
       MatTrial = cell2matrix( Y(1:end,trial)); % Cell2Matrix substitutes empty rows in the cell-array
                                                % with NaN rows
       NotNaN = find( sum(isnan(MatTrial.')) == 0); % It finds the positions of empty ROWS
       % Interpolation only on channels without NaNs
       % Empty cells [] also for decimated cell structure for not used channels
       
       if ( trial > par.MaxPlots)
           par.PlotSignalsON = 0;
       end
       
       X  = decimate_and_filtering_Matlab( MatTrial(NotNaN,:), Fs, Decim, par );
       Xcell(NotNaN,trial) = matrix2cellarray(X);
    end
    clear  Y
    eval( [ window_name{win},'= Xcell;']);
end

fs_original = fs;
Fdec        = fs/Decim;
fs          = Fdec;

%% Generate .mat file
date_time = datestr(datetime, 'mm-dd-yyyy HH-MM');
user_name = getenv('username');
% filename_out = char([filename, '-', type, '-', 'generated-on-', date_time, '-by-', user_name]);
if par.AC_removal == 1; % 1 to Remove AC by Notch-filtering at 60 Hz
 str_AC = 'ACnotch';
else
 str_AC = '';
end
if par.LP_Filter  == 1; % 1 to Activate Low-Pass filter with MANDATORY
str_LP = ['LP',num2str(par.LP_Freqs(1))];
else
 str_LP = '';
end
   
if par.BP_Filter  == 1; % for Active BP filter
str_BP = ['BP',num2str(par.BP_Freqs(2)),'-',num2str(par.BP_Freqs(3))];
else
 str_BP = '';
end

if par.HP_Filter  == 1; % 1 to Activate Low-Pass filter with MANDATORY
str_HP = ['HP',num2str(par.HP_Freqs(2))];
else
 str_HP = '';
end

str_filt = [str_AC, str_LP, str_BP, str_HP];

filename_out = extractBefore(filenameIN,type_trial);
Decim_eff = round(40000/fs); % decimation wrt to original data at 40000 Hz
filename_out = char([filename_out, 'AfterSplit-Dec', int2str(Decim_eff), str_filt, 'Nchan',int2str(Nch), '-', type_trial, '-', date_time,'.mat' ]);

clear X Xcell MatTrial Nch Ntrial NotNaN par.saveON window_name type_trial 
clear str_AC str_LP str_BP str_HP str_filt;  

save(filename_out, '-v7.3')

