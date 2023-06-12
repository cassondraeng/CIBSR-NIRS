%% Script to analyze GNG -- CIBSR fNIRS Script
% What you will need before run this script: NIRS Analyzer Toolbox, QTNirs
% toobox, and customized functions
%
% Adapted by Dr. Rihui Li & Dr. Cassie Eng
% 2-1-2023


clc; clear all; close all

%% load data
dataDir = '/Users/cassondraeng/Documents/MATLAB/BREXFIT/Rest';
raw = nirs.io.loadDirectory(dataDir);


%% ----------- Start the preprocessing -----------%%

%% Draw the raw data 
raw(1).draw


%% Downsample data (optional)
job = nirs.modules.Resample;
job.Fs = 2;  % Sets the new sample rate to 2Hz
raw = job.run(raw);

%% First, label the short channels
job = nirs.modules.LabelShortSeperation();
job.max_distance = 15;  % distance less than 15mm will be labeled as short channels
raw = job.run(raw);

%% OD conversion
j=nirs.modules.OpticalDensity; % Data conversion job declaration
od=j.run(raw);

%% Motion artifact correction using Spline interpolation and Wavelet
% jSP=nirs.modules.Homer_SplineFilter(); % job declaration
% odSP=jSP.run(od); %  job implementation
% % 
% jHWL=nirs.modules.WaveletFilter(); % job declaration
% odHWL=jHWL.run(odSP); %  job implementation


% %% Motion artifact correction using TDDR

jHWL=nirs.modules.TDDR(); % job declaration
odHWL=jHWL.run(od); %  job implementation



%% Bandpass filter
job = eeg.modules.BandPassFilter();
job.lowpass = 0.5; % use low-pass to remove cardiac and instrumental noise
job.highpass = 0.05; % Use high-pass filter to remove drift in signal
job.do_downsample = false;
odBP = job.run(odHWL);



%% Compute the HB data
job = nirs.modules.BeerLambertLaw;
job = nirs.modules.ExportData(job);
job = nirs.modules.KeepTypes(job); % Remove deoxy-hemoglobin
    job.types = {'hbo'}
hb = job.run(odBP);


%% Get the short channels
ssChan = find(hb(1).probe.link.ShortSeperation == 1);

% set all bad short channels to 0, GLM can not process short channel with nan values
cleanhb = hb;
sChannel = ssChan;
for ii = 1:numel(cleanhb)    
    cleanhb(ii).data(:,sChannel(isnan(cleanhb(ii).data(1,sChannel)))) = 0;
end


%% Modify stim information

newhb = cleanhb;

% We first remove un-related conditions
%channel 4=intertrial crosshairs; channel 8=start; channel16=end
job=nirs.modules.DiscardStims;
job.listOfStims={'stim_channel4';'stim_channel8'; 'stim_channel16'};
newhb=job.run(newhb);

% And rename the remaining stims 
% channel 2 = go trials; channel 1 = no/go trials;
job=nirs.modules.RenameStims;
job.listOfChanges={'stim_channel1'  'Go';...
    'stim_channel2' 'NoGo'};
newhb=job.run(newhb);



%% ---------- Run first and second level GLM --------------%%
% Note: for AR-IRLS GLM, don't run any preprocessing step (motion artifact correction, bandpass filtering)

%First level Stats
job = nirs.modules.GLM;
job.type='OLS'; 
AddShortSepRegressors = true;  % use short channels THIS IS THE RIGHT CODE**
Stats4Single = job.run(newhb);

% Draw the t-map
Stats4Single(1).draw('tstat', [-6 6], 'p < 0.05');


% Second level
j = nirs.modules.MixedEffects( );
% <http://www.mathworks.com/help/stats/wilkinson-notation.html> for more examples.
j.formula = 'beta ~ -1 + cond';  % only fixed effect here
j.dummyCoding = 'full';


%%%%THIS IS WHERE THE CODE BREAKS IF YOUR MONTAGES ARE DIFFERENT 
GroupStat = j.run(Stats4Single);

% Result visualization
GroupStat.draw('tstat', [-2 2], 'q < 0.05');
% disp(GroupStat.conditions);

% Result visualization
GroupStat.draw('beta', [-2 2]);
disp(GroupStat.conditions);



GroupStat.table()
