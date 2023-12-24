function [ Data, par ]  = decimate_and_filtering_Matlab( Y, Fs, Decim, par )
%  Data = decimate_and_filtering_Matlab( Y, Fs, Decim, par );
%
%  Function to decimate Signals Y by decimate function in Matlab
%  The 'decimate' function in Matlab has been set to work with an
%  antialiasing FIR filter with  ORDER = 128
%
%  If Decim = 1 the signals are simply filtered
%
%  The function, BEFORE DECIMATION, can also:
%  - Notch-filter out the AC at 60 Hz
%  - Low-pass filtering the signals
%  - Band-pass filtering the signals
%
%  All the filters are designd by Matlab "filterdesign"
%  (Signal Processing Toolbox) and specific functions are
%  embedded at the end of this function script.
%
%  Possibility to activate plots of signals and spectrum
%
% INPUTS
% Y:     Data matrix with signals by ROWS !!
% Fs:    Sampling Frequency
% Decim: Decimation Factor (Integer)
%
% OUTPUTS
% Data:     Decimated Data matrix with signals by ROWS (!)
% par       parameters at the Input = some extra string parameters to save the file
%
% OPTIONAL parameters-------
% [par.AC_removal]: [0]/1  remove AC by Notch-filtering at 60 Hz
%
% [par.LP_Filter]:  [0]/1  Activate Low-Pass filter with MANDATORY
%  par.LP_Freqs:  [Fpass, Fstop] frequencies for the filter
%
% [par.BP_Filter]:  [0]/1  Activate Band-Pass filter with MANDATORY
%  par.BP_Freqs:  [Fstop1 , Fpass1 , Fpass2, Fstop2]
%
% [par.HP_Filter]:  [0]/1  Activate High-Pass filter with MANDATORY
%  par.HP_Freqs:  [Fstop , Fpass] frequencies for the filter
%
% [par.PlotSignalsON]:   [0]/1  Plot original and decimated signals for comparison
% [par.PlotSignalsIndex]:   [1,2]  Index numbers of the Signals be plot on a 'subplot' figure
%                        (At most it will plot the first 10);
%
% [par.MaxSamplesPlot]: [inf]  Max # of original time-samples to be plot (inf = ALL)
% [par.MaxFreqPlot]:    [Fs/(2*Decim)] Max frequency shown in the plot
%
% [par.saveON]:   [0]/1 to save the decimated matrix to a file
% [par.filename]: ['DecimTime'] starting part of the name of the output file

if nargin < 4
    par.AC_removal     = 0;
end
if ~isfield(par,'AC_removal')
    par.AC_removal=0;
end
if ~isfield(par,'f_AC')
    par.f_AC = 60;
end
% if ~isfield(par,'Notch_BW')
% par.Notch_BW  = 1;
% end
if ~isfield(par,'BP_Filter')
    par.BP_Filter = 0;
end
if ~isfield(par,'LP_Filter')
    par.LP_Filter = 0;
end
if ~isfield(par,'HP_Filter')
    par.HP_Filter = 0;
end
if ~isfield(par,'saveON')
    par.saveON = 0;
end
if ~isfield(par,'filename')
    par.filename  = 'DecimTime';
end

if ~isfield(par,'PlotSignalsON')
    par.PlotSignalsON   = 0;
end
if ~isfield(par,'PlotSignalsIndex')
    par.PlotSignalsIndex =[1,2];
end

if ~isfield(par,'MaxPlots')
    par.MaxPlots = 2;
end
if ~isfield(par,'MaxSamplesPlot')
    par.MaxSamplesPlot = inf;
end
if ~isfield(par,'MaxFreqPlot')
    par.MaxFreqPlot    = Fs/(2*Decim);
end
if par.BP_Filter == 1
    if ~isfield(par,'BP_Freqs')
        error('[Fstop1,Fpass1,Fpass2,Fstop2] must be set in ''par.BP_Freqs''')
    end
end
if par.LP_Filter == 1
    if ~isfield(par,'LP_Freqs')
        error('[Fpass,Fstop] must be set in ''par.LP_Freqs''')
    end
end
if par.HP_Filter == 1
    if ~isfield(par,'HP_Freqs')
        error('[Fstop,Fpass] must be set in ''par.HP_Freqs''')
    end
end

MaxPlots     = 10;
[NumChan, N] = size(Y); % Data by Rows
% If Y is a vector, it forces Y to be a row vector
if min([NumChan, N])==1
    if  N == 1
        Y = Y.';
        [NumChan, N] = size(Y);
    end
end

NumSubPlots = min([MaxPlots,NumChan,length(par.PlotSignalsIndex)]);

if par.AC_removal == 1
    [b_AC,a_AC] = Notch_filter_design(Fs);
    Xfil = filtfilt(b_AC, a_AC, Y.');
    str_AC ='AC-notch';
else
    Xfil = Y.';
    str_AC ='';
end

if par.BP_Filter == 1
    [b_BP,a_BP] = BP_filter_design(Fs, par);
    Xfil = filtfilt(b_BP, a_BP, Xfil);
    
    F1 = int2str(par.BP_Freqs(2));
    F2 = int2str(par.BP_Freqs(3));
    str_BP =[ '-BP', F1,'_' F2 ];
else
    str_BP = '';
end

if par.LP_Filter == 1
    [b_LP,a_LP] = LP_filter_design(Fs, par );
    Xfil = filtfilt(b_LP, a_LP, Xfil);
    F1     = int2str(par.LP_Freqs(1));
    str_LP =[ '-LP', F1];
else
    str_LP = '';
end

if par.HP_Filter == 1
    [b_HP,a_HP] = HP_filter_design(Fs, par );
    Xfil = filtfilt(b_HP, a_HP, Xfil);
    F1     = int2str(par.HP_Freqs(2));
    str_HP =[ '-HP', F1];
else
    str_HP = '';
end


L = ceil(N/Decim);    % # Samples of the decimated signal
Data = zeros(NumChan,L); % Preallocation of the output DataMatrix
% Downsampled signal in time-domain

if Decim == 1
    Data = Xfil.';
else
    for ind =1:NumChan
        if sum(isnan(Xfil(:,ind)))> 0
            Data(ind,:)= NaN*ones(1,L);
        else
            Data(ind,:) = decimate( Xfil(:,ind), Decim, 128, 'FIR'); % 128th order FIR anti-aliasing filter
        end
    end
end

% plot_num = 0;
if (par.PlotSignalsON == 1) % && (plot_num <= par.MaxPlots)
    % plot_num = plot_num+1;
    for ind =1:NumSubPlots
        if ind == 1
            fig1 = figure;
            fig2 = figure;
        end
        % t     = (0:Npoint-1)/Fs;       % Original time axis
        % t_dec = (0:Decim:Npoint-1)/Fs; % Downsampled time axis
        % f     = (0:M-1)*DeltaF;        % frequency axis
        Nplot  =  min([par.MaxSamplesPlot, N]); % Number of points plotted in time
        
        % Time domain plots
        figure(fig1); subplot( NumSubPlots, 1, ind);
        
        plot( (0:Nplot-1)/Fs, Y(ind,1:Nplot),'k'); hold on;
        plot( (0:Decim:Nplot-1)/Fs, Data(ind,1:length((0:Decim:Nplot-1))),'r'); grid on;
        xlabel('time (sec.)'); ylabel('Time');
        legend( 'Original', 'Decimated');
        title(['Channel #', int2str( par.PlotSignalsIndex(ind) )]);
        
        
        % if mod(L,2)~=0
        %     warning('FFT_length / Decimation is not an even integer number');
        %     warning('FFT_length increased to have integer even L')
        %     L_2 = ceil(0.5*N/Decim); % Integer half-duration of downsampled signal
        %     L   = 2*L_2;             % Even duration of dowsampled signal
        %     M   = L*Decim;           % Zero-padded signal length (DFT size)
        % else
        %     L_2 = 0.5*L;  % Integer Half duration
        %     M    = N;     % Original signal length
        % end
        
        %     L_2 = 0.5*L;  % Integer Half duration
        yf     = fft(Y(ind,:))/Fs;
        DeltaF = Fs/N;
        N_2    = floor(0.5*N);
        xf     = fft(Data(ind,:))/(Fs/Decim);
        DeltaF2 = Fs/(Decim*L);
        L_2     = floor(0.5*L);
        
        % Frequency domain plots
        figure(fig2); subplot(NumSubPlots,1, ind);
        plot( (0:N_2-1)*DeltaF, 20*log10(abs(yf(1:N_2))),'k'); hold on;
        plot( (0:L_2-1)*DeltaF2, 20*log10(abs(xf(1:L_2))),'r'); grid on;
        AX = axis;
        axis([0,par.MaxFreqPlot, AX(3),AX(4)]);
        
        xlabel('frequency (Hz)'); ylabel('|X(f)|^2 (dB)')
        legend( 'Original', 'Decimated');
        title(['Channel #', int2str(par.PlotSignalsIndex(ind))]);
    end
end

Fdec    = Fs/Decim;
Duration_sec = N/Fdec; % Duration in seconds
str_sec      = num2str(Duration_sec,'%.1f');
str_dec      = int2str(Decim);

% Save in a.mat file the data matrix with signals per rows.
% The file name encodes some parameters
% (Decimation, AC-removal, Duration in seconds of the series)

str_filt = [ str_AC, str_BP, str_LP , str_HP ];
par.filename = [par.filename ,str_filt,'-Dec', str_dec, '-Nchan' int2str(NumChan), '-Seconds' str_sec '.mat' ];

parfilters = par.LP_Filter + par.HP_Filter + par.BP_Filter + par.AC_removal;
filters_ON = (parfilters ~= 0);
decim_ON   = (Decim ~= 1);
save_ON    = (par.saveON ==1);
if ( save_ON && (filters_ON || decim_ON ))
    % -v7.3 option can be removed if the file size is not huge
    eval( ['save ', par.filename , ' Data Fdec Fs par -v7.3'] );
end

% FUNCTIONS to design filters
% Some parameters are fixed but they can be edited or transformed as
% input parameters

function [b,a,Hd] = Notch_filter_design(Fs)
%FILTERDESIGNNOTCHCHEBYSHEV
% Returns ALSO a discrete-time filter object Hd.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and Signal Processing Toolbox 8.2.
% Generated on: 06-Apr-2020 10:53:50

% Chebyshev Type I Bandstop filter designed using FDESIGN.BANDSTOP.

% All frequency values are in Hz.
% Fs = 256;  % Sampling Frequency

Fpass1 = 59;          % First Passband Frequency
Fstop1 = 59.5;        % First Stopband Frequency
Fstop2 = 60.5;        % Second Stopband Frequency
Fpass2 = 61;          % Second Passband Frequency
Apass1 = 0.5;         % First Passband Ripple (dB)
Astop  = 30;          % Stopband Attenuation (dB)
Apass2 = 0.5;         % Second Passband Ripple (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.bandstop(Fpass1, Fstop1, Fstop2, Fpass2, Apass1, Astop, ...
    Apass2, Fs);
Hd = design(h, 'cheby1', 'MatchExactly', match);

% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);

% [EOF]

function [b,a,Hd] = BP_filter_design(Fs, par)
%FILTERDESIGNBPELLIPTIC70_110
% Returns ALSO a discrete-time filter object Hd.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and Signal Processing Toolbox 8.2.
% Generated on: 06-Apr-2020 11:01:13

% Elliptic Bandpass filter designed using FDESIGN.BANDPASS.
% All frequency values are in Hz.
% Fs = 256;  % Sampling Frequency

% Fstop1 = 68;  % First Stopband Frequency
Fstop1 = par.BP_Freqs(1);
% Fpass1 = 70;  % First Passband Frequency
Fpass1 = par.BP_Freqs(2);
% Fpass2 = 110;     % Second Passband Frequency
Fpass2 = par.BP_Freqs(3);
% Fstop2 = 112;     % Second Stopband Frequency
Fstop2 = par.BP_Freqs(4);

Astop1 = 40;      % First Stopband Attenuation (dB)
Apass  = 0.25;    % Passband Ripple (dB)
Astop2 = 40;      % Second Stopband Attenuation (dB)
match  = 'both';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
    Astop2, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);

% [EOF]

function [b,a,Hd] = LP_filter_design(Fs, par)
%FILTERDESIGNLOWPASSCHEBYSHEV70
% NOTE: by using "filtfilt" the atteniations are DOUBLED in dB
% Returns ALSO a discrete-time filter object Hd.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and Signal Processing Toolbox 8.2.
% Generated on: 06-Apr-2020 11:24:00

% Chebyshev Type I Lowpass filter designed using FDESIGN.LOWPASS.

% All frequency values are in Hz.
% Fs = 256;  % Sampling Frequency

% Fpass = 68;
Fpass = par.LP_Freqs(1); % Passband Frequency
% Fstop = 70;
Fstop = par.LP_Freqs(2); % Stopband Frequency
Apass = 0.25;         % Passband Ripple (dB)
Astop = 40;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its CHEBY1 method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
Hd = design(h, 'cheby1', 'MatchExactly', match);

% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);

% [EOF]

function [b,a,Hd] = HP_filter_design(Fs, par)
%FILTERDESIGNHIGHPASSELLIPTIC0.5HZ Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.6 and Signal Processing Toolbox 8.2.
% Generated on: 10-Apr-2020 11:01:28

% Elliptic Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
% Fs = 256;  % Sampling Frequency

% Fstop = 0.4;       % Stopband Frequency
Fstop = par.HP_Freqs(1);   % Stopband Frequency
% Fpass = 0.5;        % Passband Frequency
Fpass = par.HP_Freqs(2);   % Stopband Frequency
Astop = 20;          % Stopband Attenuation (dB)
Apass = 0.05;        % Passband Ripple (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);



function Hd = HPfilterNEW1000Hz
%HPFILTERNEW1000HZ Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.8 and Signal Processing Toolbox 8.4.
% Generated on: 01-Jul-2020 20:07:43

% Elliptic Highpass filter designed using FDESIGN.HIGHPASS.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fstop = 0.3;     % Stopband Frequency
Fpass = 0.5;     % Passband Frequency
Astop = 40;      % Stopband Attenuation (dB)
Apass = 0.1;     % Passband Ripple (dB)
match = 'both';  % Band to match exactly

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
Hd = design(h, 'ellip', 'MatchExactly', match);

% Get the transfer function values.
[b, a] = tf(Hd);

% Convert to a singleton filter.
Hd = dfilt.df2(b, a);



% [EOF]



% [EOF]


