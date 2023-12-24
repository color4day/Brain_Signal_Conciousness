%%%

Decim = 1; % Decimation Factor (1  if already decimated)

% Filtering options
PAR.AC_removal = 0;   % 1 to Remove AC by Notch-filtering at 50 Hz

PAR.LP_Filter  = 1;   %   1 to Activate Low-Pass filter with MANDATORY
PAR.LP_Freqs   = [30,35]; % [Fpass, Fstop] frequencies for the filter
  
PAR.BP_Filter  = 0;       % 1 to Activate Band-Pass filter with MANDATORY
PAR.BP_Freqs   = [14,15,30,31]; % [Fstop1 , Fpass1 , Fpass2, Fstop2]

PAR.HP_Filter = 0;         % 1 to Activate High-Pass filter with MANDATORY
PAR.HP_Freqs  = [0.3,0.5]; % [Fstop,Fpass] frequencies for the filter
                           % Elliptic Filter designed specifically to High-Pass
                           % 0.5 Hz, Fs = 256Hz, 40dB stopband with "filtfilt"
                           % See details in the function 
PAR.PlotSignalsON = 1;     % 1 to plot original and decimated and filtered signals for comparison
PAR.MaxPlots      = 2;     % maximum number of Trials to be plotted
PAR.PlotSignalsIndex = [1,21]; %,33,47,72,86,97,111]; % Index numbers of signals and spectrum to be plotted
                           % on a 'subplot' figure
                           % (Only the first 10 will be actually plotted);

PAR.MaxSamplesPlot = 300 ;      % Max # of original time-samples to be plotted (inf = ALL)
% par.MaxFreqPlot = Fs/(2*Decim); % Max frequency shown in the plot 

DirStruct = dir('*.mat');
% FileNames = char( {DirStruct.name}.' );

% filenameIN ='M_20180530-orientation-north-generated-on-05-04-2020 16-48.mat';
% type_trial = 'orientation-north'; % which MUST be contained in finenameIN;
% funcDecimateMonkeyDataset_after_Splitting( filenameIN, type_trial, Decim);

for ind = 1:1 %size(DirStruct,1)
    filename = DirStruct(ind).name;
    VAR = load(filename); % load all the variables stored 
    PAR.MaxFreqPlot = VAR.fs/(2*Decim); % Max frequency shown in the plot (Check fs !!)
    funcDecimateMonkeyDataset_after_Splitting2_woChoice( filename , VAR.type, Decim, PAR);
end 
