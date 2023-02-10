%% Script to analyze brexfit data--CIBSR
clc; clear all; close all

% This script does no one analysis in particular.
% It is the first script I started writing and has turned into a 
% collection of notes, resources, code, and different approaches that might
% come handy. 


%% recursive search subdirectories for 
% see this link for explanation
% https://se.mathworks.com/matlabcentral/answers/429891-how-to-recursively-go-through-all-directories-and-sub-directories-and-process-files
% possibly see this (to make it load properly) https://blogs.mathworks.com/pick/2016/12/16/glob-file-searching-in-matlab/

rootdir = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS';
%filelist = dir(fullfile(rootdir, '**\*.*'))  %get list of files and folders in any subfolder
nback_dirlist = dir(fullfile(rootdir, '**\*_nback\2022*'));
%3\Data\NIRS\3004\indoor\in_post_gonogo\2022-07-18_006\'
nback_dirlist.folder; %lists all folders 

raw = [] % declare empty data set 
for i=1:numel(nback_dirlist)
    filedir = nback_dirlist(1).folder; % find file directory
    new_file = nirs.io.loadDirectory(filedir); % load file
    
    % check that files are nback, filter files here if not
    name = new_file.demographics.experiment;
    startidx = strfind(name, '_nback');
    if isempty(startidx)
        disp('THis is not nback!/ have a look at this file:')
        disp(new_file.description)
    end

    % Possibly put this in an elseif to add only files that "pass" the "is
    % this nback?"-test
    raw = [raw new_file] % append new files to data set
 
end

raw.draw()

% stimuli does not appear to be properly loaded this way
% could load .nirs file spesifically
%%
% Note: It is MUCH more efficent to use loadDirectory than to use loadSNIRF 
% (using my computer (HD) and the attempt below). 

% qr_dirlist = dir(fullfile(rootdir, '**\qr_during\2022*\*.snirf')) % try to get snirf files
% qr_dirlist.name % gives the name of the file
% qr_dirlist.folder % gives the path to the folder
% % --> must combine the two above to get the full path to the file
% % Either way, there must combine two strings with this approach
% raw_qr = [] % declare empty data set
% for i=1:1:length(qr_dirlist) % iterate through dirlist 
%     path = append(qr_dirlist(i).folder, '\', qr_dirlist(i).name) % find full path of individual sfnirs files
%     raw = nirs.io.loadSNIRF(path) % load individual sfnirs files
%     raw_qr = [raw_qr raw] % append individual files to data set 
% end

%% load data

% try loading whole directory. that works, loads alll files
datadirall = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS'
datadir1 = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS\3004'
datadir2 = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS\3023'

raw = nirs.io.loadDirectory(datadir1)

% potential for error!! different montages are used for cycling exercise,
% these must be re-registered to another probe o loaded separately 
% think this depends also on the order of which files are loaded, cycling
% vs not

% Loading NIRx file geometry from:
%      Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS\3004\indoor\in_post_gonogo\2022-07-18_006\
%       Note: This registration will be used for all subjects

% visual inspection 

raw.draw()

raw(1).description
figure; raw(1).draw()
raw(1).stimulus

raw(2).description
figure; raw(2).draw()
raw(2).stimulus

raw(3).description
figure; raw(3).draw()
raw(3).stimulus

raw(4).description
raw(4).demographics.experiment

figure; raw(3).draw()
raw(3).stimulus

demoTable = nirs.createDemographicsTable(raw)
% Save table!

%% split data set according to task
% need to separate n-back, gonongo and during. either by loading the
% correct files only or pruning the dataset after

% start with nback

raw_nback = [] % declare empty dataset for nback data only

for i=1:numel(raw)
    name = raw(i).demographics.experiment;
    startidx = strfind(name, '_nback');
    if isempty(startidx)
        disp('THis is not nback!')
    elseif isempty(startidx) == 0
        disp('This is nback')
        % append to nback data set
        raw_nback = [raw_nback raw(i)]
    end
end

% visual inspection of nback data
raw_nback.draw()

% some are missing stimuli markers and i dont know why 
raw_nback(1).demographics.experiment
figure; raw_nback(1).draw()

raw_nback(6).demographics.experiment
figure; raw_nback(6).draw()
raw_nback(6).description

% the other have stimuli markers
raw_nback(3).demographics.experiment
figure; raw_nback(3).draw()
raw_nback(3).description

raw = raw_nback  


%% plot montage

figure; raw_nback(1).probe.draw()
raw_nback(1).probe.defaultdrawfcn='10-20'
raw_nback(1).probe.defaultdrawfcn='?'
raw_nback(1).probe.defaultdrawfcn='3D mesh'

%% correct stimuli name and duration 

% discard stimuli not in use
job=nirs.modules.DiscardStims;
job.listOfStims={
    'stim_aux1'; 
    'stim_aux2'; 
    'stim_aux3'; 
    'stim_aux4'};
raw=job.run(raw);

% check that it works
figure; raw(2).draw()

% rename stimuli --> Uncomment and replace with correct names once
% ascertained
% j = nirs.modules.RenameStims()
% j.listOfChanges = {
%     'stim_channel1', 'name?';
%     'stim_channel2', 'name';
%     'stim_channel4', 'namememhehhehehe';
%     'stim_channel8', 'anothernameme';
%     'stim_channel16', 'wooowolongmname'
%     };
% raw = j.run(raw);

% Looks like 8 and 16 are the ones with actual stimuli here

raw(3).stimulus

T = nirs.createStimulusTable(raw)

%% ?? change stimuli duration???? 
% change stimulus durations
%rawChanged = nirs.design.change_stimulus_duration(raw,{'stim_channel1', 'stim_channel2','stim_channel3'},4);
rawChanged = nirs.design.change_stimulus_duration(raw(3),{'stim_channel1', 'stim_channel2','stim_channel4','stim_channel8','stim_channel16' },4);
figure; rawChanged.draw()


%% label the short-channels
job=nirs.modules.LabelShortSeperation();
job.max_distance=15;
raw = job.run(raw);

% check that short-channels are labeld (command might produce error if the
% field ShortSeparation does not exist
raw(1).probe.link.ShortSeperation == 1 % some 1s should appear

%% quality check of raw data  - QUALITY assessment via QTNIRS job
% is this optional ?????
% the results from this quality check is not really used for anything right
% now (but it might be later), code works and it might come handy after
% having read the paper associated with QTNIRS

j = nirs.modules.QT()
j.qThreshold =0.75
j.fCut = [0.5 2.0]
ScansQuality = j.run(raw)
for ii = 1:numel(ScansQuality)
   figure(ii); ScansQuality(ii).drawGroup('sq');
    
end

% Look at the data
ScansQuality(1).qMats
ScansQuality(1).qMats.bad_links
ScansQuality(1).qMats.good_combo_link
ScansQuality(1).probe.link


%% preproccesing: conversion to hbo/hbr

% basic preproccesing
j = nirs.modules.default_modules.basic_preprocessing();
list = nirs.modules.pipelineToList(j);
disp(list);
%  {1×1 nirs.modules.ImportData    }
%     {1×1 nirs.modules.RemoveStimless}
%     {1×1 nirs.modules.FixNaNs       }
%     {1×1 nirs.modules.Resample      }
%     {1×1 nirs.modules.OpticalDensity}
%     {1×1 nirs.modules.BeerLambertLaw}
%     {1×1 nirs.modules.TrimBaseline  }
%     {1×1 nirs.modules.ExportData    }
list{4}.Fs = 1;
list{7}; % Trim baseline, leave defaults of 30 s pre/post stimuli
j = nirs.modules.listToPipeline(list) % convert back to pipeline
hb = j.run(raw);

hb.draw() % visual inspection

%% 1-st level stats

% option 1
% motion artifact correction by using AR-IRLS.
j = nirs.modules.GLM();
SubjStats = j.run(hb);

% option 2
% short-channel regressors,
% motion artifact correction by using AR-IRLS.
j = nirs.modules.GLM();
j.AddShortSepRegressors = true;
SubjStatsSS = j.run(hb);

% option 3
% short-channel regressors,
% motion artifact correction by including accelerometer as a regressor 
% motion artifact correction by using AR-IRLS.
% NOTE: Must inspect accelerometer data visually and remove the
% accelerometer that (most often) isn't placed on participants head, and
% thus isn't giving accurate data
j = nirs.modules.AddAuxRegressors();
j = nirs.modules.GLM(j);
j.AddShortSepRegressors = true;
SubjStatsSSaux = j.run(hb);

%% compare the results of the various options 
% I expect the betas to differ
% currently the don't, but this is only one participant and the data is
% resampled at a low frequency
head(SubjStats(1).table())

% draw results 
SubjStats(1).table('q <0.05')

figure; SubjStats(1).draw('tstat', [-5 5], 'q <0.05')

figure; SubjStats(1).draw('tstat', [-5 5], 'q <0.05')


% search for values in table
for i=1:numel(SubjStats)
    disp('=============================================================')
    StatsTable = SubjStats(i).table();
    head(StatsTable)
    StatsSSTable = SubjStatsSS(i).table();
    head(StatsSSTable)
    StatsSSauxTable = SubjStatsSSaux(i).table();
    head(StatsSSauxTable)
    disp('=============================================================')
%     idx = find((StatsTable.q <0.05));
%     %idx = find((StatsTable.q <0.05) & strcmp(StatsTable.cond,'RampGame3'));
%     disp(StatsTable(idx, :)); 
end

%% subject/participant leverage



%% 2-nd level analysis

j = nirs.modules.MixedEffects();
% We must specify the formula for the mixed effects.  This one calculates
% the group mean for each condition.  There is also a random intercept for
% each subject.  Google "matlab wilkinson notation" or see
% <http://www.mathworks.com/help/stats/wilkinson-notation.html> for more examples.
j.formula = 'beta ~ -1 + cond + (1|subject)';
j.dummyCoding = 'full';
GroupStats = j.run(SubjStats);

% (in final analysis) change defaults to: run robust stats, include diagnostics, and print
% what's happening while the model runs
% j.robust = true;
% j.include_diagnostics=true;
% j.verbose = true;


% Conditions relative to baseline (Beer-Lambert law baseline)
% Draw the probe using t-stat & false discovery rate (q value)
GroupStats.draw('tstat', [], 'q < 0.05')

%% contrasts
disp(GroupStats.conditions);

% define contrast(s) 
c = {'stim_channel8-stim_channel16'}

% Calculate stats with the ttest function
ContrastStats = GroupStats.ttest(c);

head(ContrastStats.table)
idx = find(ContrastStats.q < 0.05); 
SignificantContrastStats = ContrastStats.table; 
SignificantContrastStats = SignificantContrastStats(idx, :); 
disp(SignificantContrastStats);
%ContrastStats.draw('tstat',[-10 10],'q<0.05')
ContrastStats.draw('tstat',[],'q<0.05')


%% HRF 

% HRF?
HRF = GroupStats.HRF
nirs.viz.plot2D(HRF)
HRF.draw()

%% suggestions for next steps?

% import files from more than one participant
% import files using a glob/recurrsive search 
% import .nirs files (perhaps in combination with the one above)
% change stimuli duration (this is important, need to work, look ath the
% code form the load paper you already have maudenbauer something, and lok
% at hrf nback paprs and papers citing the nirs-toolbox

% Anticipation of a mentally effortful task recruits Dorsolateral
% Prefrontal Cortex: An fNIRS validation study
% 10.1016/j.neuropsychologia.2018.04.033











