%% Script that analyzes gonogo fnirs data from Brexfit
% To run this script you need nirs-toolbox and QT-nirs 

% Last updated: 2023-01-23
% Created by: Henrikke Dybvik
% Using MATLAB R2021b

% Question: Will we use iti and rest markers?
% TODO: Decide on group level model + contrasts
% HD: I want to explore the hrf response options. 
% Add piece of code that removes channels based on QC-NIRS results.

% Flag channels with less than 10% high qulity channels? 

%% load gonogo data 
% Uncomment/comment dirlist and rootdir below according to os:

% [For Windows]
rootdir = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS';
dirlist = dir(fullfile(rootdir, '**\*_gonogo\202*')); % get folders with gonogo data

% [For Mac]
%The rootdir needs to match how you login to the server, which is currently
%via Cisco Anyconnect with the following path entered into the "Go"
%connection smb://dibshome/Projects/NIRS_Projects/SpecializedFoundation/BREXFIT_S3
%rootdir = '/Volumes/BREXFIT_S3/Data/NIRS/';
%dirlist = dir(fullfile(rootdir, '*/*/*_gonogo/202*')); % get folders with gonogo data
%dirlist.name
%dirlist.folder

raw = []; % declare empty data set
for i=1:1:length(dirlist) % iterate through dirlist 
    path = append(dirlist(i).folder, '\', dirlist(i).name); % [Windows] find full path of individual folders
    %path = append(dirlist(i).folder, '/', dirlist(i).name); % [Mac] find full path of individual folders
    raw_indi = nirs.io.loadDirectory(path); % load individual folders/files
    raw = [raw raw_indi]; % append individual files to data set 
end

rawgonogo = raw
%% Search for a particular file or participant
for i =1:length(raw)
    if strcmp(raw(i).demographics.subject, '3012')
        raw(i)
        figure; raw(i).draw()
    end
end

for i =1:length(hb)
    if strcmp(hb(i).demographics.subject, '3012')
        hb(i)
        figure; hb(i).draw()
    end
end

%% Filter data based on QC nirs 

% Comments from NIRS_QC excel sheet: 
% Participant 3004 for in_post_gognogo: triggers were not appearing for 
% gonogo…nirs 6.2 mb. --> HD: I don't think there is much we can do about
% that. Thus, no additional filtering/data selection needed at this moment
% (written 2022-11-29). 

demoTablegonogo = nirs.createDemographicsTable(rawgonogo);
disp(demoTablegonogo)

%% Add demographics filed according to condition 
% We have the variable timing (pre, post) and the variable load (qr/out/in)
% These can be placed in demographics table and be used as predictors or 
% covariates in the group level model later. 

for i =1:numel(raw)
    % rename experiment variable 
    newExpStr = erase(raw(i).demographics.experiment, '_gonogo');
    newExpStr = strip(newExpStr,'left','_'); % remove leading underscores from string
    raw(i).demographics.experiment = newExpStr;
    
    % create timing variable (indicates pre/post)
    if endsWith(newExpStr,'post') | endsWith(newExpStr,'postrs') | endsWith(newExpStr,'post_restart')
        raw(i).demographics.timing = 'post';
    elseif endsWith(newExpStr, 'pre')
        raw(i).demographics.timing = 'pre';

    else
        raw(i).demographics.timing = 'I donnnooooooooooooo';
    end
    
    % create load variable (indicated qr/out/in)
    if startsWith(newExpStr,'qr')
        raw(i).demographics.load = 'qr';
    elseif startsWith(newExpStr, 'in')
        raw(i).demographics.load = 'in';
    elseif startsWith(newExpStr, 'out')
        raw(i).demographics.load = 'out';
    else
        raw(i).demographics.load = 'I donnnooooooooooooo';
    end

    % convert subject id to numeric and create a character version of
    % subejct id
    raw(i).demographics.subject_id = raw(i).demographics.subject;
    raw(i).demographics.subject = str2num(raw(i).demographics.subject);
end
% create new demographics table and check that variables have been
% correctly named

demoTable = nirs.createDemographicsTable(raw);
disp(demoTable)
%% Adding demographic information
% e.g., age and baseline executive function 
j = nirs.modules.AddDemographics();
% This job contains a field (demoTable) that is used to populate the
% demographics.  This can be programatically or simply read from a CSV file

additional_demographics = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\Demographics\BrainResearchExercis-Demographics2.xlsx';
j.demoTable = readtable(additional_demographics);% [root_dir filesep 'demo_data' filesep 'data' filesep 'demographics.csv'] );

% We are going to match the subject column in the above table
j.varToMatch = 'subject';
% Note "subject" (or whatever you are matching based on) needs to be an
% entry in BOTH the CSV file and the nirs.core.<> class that you are
% loading the demographics to.
% subject must be numeric for this to work

raw = j.run(raw); %SubjStats);  Can run it on subjstats as well. 

%% add age for both child and adult to the same demographic varaible 

for i =1:numel(raw)
    if (isnan(raw(i).demographics.age_child)) & (~isnan(raw(i).demographics.age_adult))
        %use adult age
        raw(i).demographics.age = raw(i).demographics.age_adult;
    elseif (~isnan(raw(i).demographics.age_child)) & (isnan(raw(i).demographics.age_adult))
        % use child age
        raw(i).demographics.age = raw(i).demographics.age_child;
    else
        disp('something is wrong');

    end
end
%% Order dataset so pre occurs first (and hence will be used as reference 
% category in group level mixed model later on)
%
% Can adjust reference group

for i=1:numel(raw)
    if strcmp(raw(i).demographics.timing,'pre') % Change the string ('text in brackets') to decide what will be used as reference in group level model 
        % identify first occurence of pre
        rawpostfirst = raw(i);
        idx_remove = i;
        break
    else
        continue
    end
end
raw(idx_remove) = [];
raw = [rawpostfirst raw];


%% Visual inspection
% raw.draw() % draw the entire data set 
figure; raw(4).draw()

demoTable_gonogo = nirs.createDemographicsTable(raw);
disp(demoTable_gonogo)

% Check which stimuli we have?
%raw.stimulus
raw(4).stimulus

%% correct stimuli name and duration 

% discard stimuli not in use
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_aux1'; 
    'stim_aux2'; 
    'stim_aux3'; 
    'stim_aux4'};
raw=j.run(raw);

% check that it works
%figure; raw(3).draw()

j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_channel4'; 
    'stim_channel8'; };
raw=j.run(raw);

% Event markers (from task script):
% TARGET_TRIAL = 1
% NON_TARGET_TRIAL = 2
% ITI_TRIAL = 4
% REST = 8
% MID_TASK_BREAK = 16   This marker does not seem to be in the data

% rename stimuli
j = nirs.modules.RenameStims()
j.listOfChanges = {
    'stim_channel1', 'noGo'; % X
    'stim_channel2', 'Go'; % every other letter
    'stim_channel4', 'iti';
    'stim_channel8', 'rest';
    };
raw = j.run(raw);

% Question: Will we use iti and rest markers? Nope, not for now. 

% discard additional stimuli if neccessary
j=nirs.modules.DiscardStims;
j.listOfStims={
    'iti'; 
    'rest'; };
raw=j.run(raw);

% remove files without simuli 
% print how many files 
fprintf('number of files in data set before removing stimless files: %i\n', numel(raw));
j = nirs.modules.RemoveStimless();
raw = j.run(raw);
fprintf('number of files in data set after removing stimless files: %i\n', numel(raw));

% lost 3 files

% label the short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15; % check units 
raw = j.run(raw);

% check that it workes
raw(1).probe.link.ShortSeperation

% Leave duration of stimuli default

%% Manual Visual inspection of 


nirs.viz.TimeSeriesViewer(hb(2))

nirs.viz.nirsviewer(hb(2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visual quality check by plotTraces
data = [];
data = hb(2).data(:,1:2:end);
%data = hb(2).data(:,1:end);
figure, hold on
for jj = 1:1:size(data,2)
   subplot(12,10,jj), plot(data(:,jj));title(num2str(jj))  
end
meanamp = mean(abs(raw(2).data(:,[1:2:end]))); figure, plot(meanamp)
figure, plot(data(:,73))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

badchannelt1{2} = [1 2 ];

% set channels to zero 
hb(1, 1).data(:,badchannelt1{2}) = []  
hb(1).probe.link(badchannelt1{2},:) = []

i=2
% only need to plot every other channel
for j = 1:2:size(raw(i).data,2)
    x = raw(i).data(:,j);
    time = raw(i).time;
    Fs = raw(i).Fs;
    
    figure; cwt(x,Fs);

    tt =  table2array(raw(2).probe.link(j,{'source', 'detector'}));
    fprintf('Figure: %i   Source: %i Detector: %i \n',j, tt(1), tt(2));
end



figure, hold on
for jj = 1:size(raw(i).data,2)
    x = raw(i).data(:,j);
    time = raw(i).time;
    Fs = raw(i).Fs;
    
   subplot(12,10,jj), 
   cwt(x,Fs);; 
   title(num2str(jj))  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hb(1).probe.link

Hb()

% set channels to zero 
hb(1, 1).data(:,1) = []  
hb(1).probe.link(1,:) = []

%% try top put it in a loop

%visual quality check by plotTraces
for i = 1:numel(hb)
    data = [];
    data = hb(i).data(:,1:2:end);
    figure, hold on
    for jj = 1:1:size(data,2)
       subplot(12,10,jj), plot(data(:,jj));title(num2str(jj))  
       %disp(figure )
    end
    meanamp = mean(abs(raw(i).data(:,[1:2:end]))); figure, plot(meanamp)
    figure, plot(data(:,73))
end

%% Quality check via QT-NIRS

j = nirs.modules.QT();
j.qThreshold =0.6;
j.sciThreshold = 0.5;
%j.pspThreshold = somevalue;
% psp default
j.fCut = [0.5 2.5];
ScansQuality = j.run(raw);

% Questions: 
% What do we choose as psp and sci threshold?

% view group results to identify loaction of bad channels on the montage. 
ScansQuality.drawGroup('sq');
ScansQuality.drawGroup('sci');
ScansQuality.drawGroup('psp');
ScansQuality.drawGroup('bar');
ScansQuality.drawGroup('sqmask');

% view individual results
for ii = 1:2 %numel(ScansQuality)
    figure;

    subplot(1,3,1)
    %figure; 
    ScansQuality(ii).drawGroup('sq')
    %title('sq')
    
    subplot(1,3,2)
    %figure; 
    ScansQuality(ii).drawGroup('sci')
    %title('sci')
    
    subplot(1,3,3) 
    ScansQuality(ii).drawGroup('psp')
   % title('psp')

    %title('Quality check:', char(raw(1,i).demographics.subject))
end

%% Identify number of bad (percentage) and their location 
% Declare variables that will hold overall survival rate
ChTotALL = 0;
badChTotALL = 0;
ChLongALL = 0;
badChLongALL = 0;
ChssALL = 0;
badChssALL = 0;

fprintf('Survival rate:\n');
for i=1:1:numel(ScansQuality)
    % find total number of channels
    Ch = raw(i).probe.link;
    ChTot = height(Ch)/2; % divide by two to account for hbo/hbr
    
    % find bad channels
    idxBadCh = find(ScansQuality(i).qMats.MeasListAct==0);
    badChTot = length(idxBadCh)/2;
    
    % find long channels
    idxChLong = find(raw(i).probe.link.ShortSeperation == 0);
    ChLong = length(idxChLong)/2;
    
    % find bad long channels
    idxBadChLong = find((ScansQuality(i).qMats.MeasListAct==0) & (raw(i).probe.link.ShortSeperation == 0));
    badChLong = length(idxBadChLong)/2;
    
    % find short channels
    idxChss = find(raw(i).probe.link.ShortSeperation == 1);
    Chss = length(idxChss)/2;
    
    % find bad short channels
    idxBadChss = find((ScansQuality(i).qMats.MeasListAct==0) & (raw(i).probe.link.ShortSeperation == 1));
    badChss = length(idxBadChss)/2;
    
    %fprintf('Scan:%i #BadChannelsTotal: %i (%.2f percent)\n',i,badChTot, (badChTot/ChTot)*100);
    %fprintf('Scan:%i #BadLongChannels: %i (%.2f percent of all long channels)\n',i,badChLong, (badChLong/ChLong)*100);
    %fprintf('Scan:%i #BadShortChannels: %i (%.2f percent of all short channels)\n',i,badChss, (badChss/Chss)*100); 
    fprintf('Scan:%i #ChannelsTotal: %.2f percent #ChannelsLong: %.2f percent #ChannelsShort: %.2f percent\n',i, (1-(badChTot/ChTot))*100, (1-(badChLong/ChLong))*100, (1-(badChss/Chss))*100);
    
    % Add individual bad channel data to overall bad channel data
    ChTotALL = ChTotALL + ChTot;
    badChTotALL = badChTotALL + badChTot;
    ChLongALL = ChLongALL + ChLong;
    badChLongALL = badChLongALL + badChLong; 
    ChssALL = ChssALL + Chss;
    badChssALL = badChssALL + badChss;
end

% Calculate survival rate overall (for all files) 
ChTotSurvALL = (1-(badChTotALL/ChTotALL))*100; % Total channels survived
ChLongSurvALL = (1-(badChLongALL/ChLongALL))*100; % Long channels survived 
ChssSurvALL = (1-(badChssALL/ChssALL))*100;

fprintf('Overall survival rate: #ChannelsTotal: %.2f percent #ChannelsLong: %.2f percent #ChannelsShort: %.2f percent\n', ChTotSurvALL, ChLongSurvALL, ChssSurvALL);


%% Preproccesing 
%Basic preproccesing pipelines
j = nirs.modules.FixNaNs();
j = nirs.modules.Resample(j);
j.Fs = 5; %resample to 5 Hz
j = nirs.modules.OpticalDensity(j); % convert from raw light intensities to optical denisty
j = nirs.modules.BeerLambertLaw(j); % convert to hemoglobin concentrations
hb = j.run(raw); % 

% visual inspection
%figure; hb(2).draw()
%figure; hb(35).draw()
%figure; hb(79).draw()

% label the short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15; % check units 
hb = j.run(hb);

%% 1-st level stats (subject level)

j = nirs.modules.Resample(j);
j.Fs = 1; %resample to 5 Hz
j = nirs.modules.GLM();
j.AddShortSepRegressors = 1; % Adding the short channel measurements as regressors to the subject level GLM (as opposed to using them as a prefilter)
%SubjStatsSub = j.run(hb(1:5)); %run on subset
SubjStats = j.run(hb);

%nirs.design.basis.Vestibular
%nirs.design.trend.constant
%j.trend_func
%j.basis
%figure; SubjStats(1).HRF.draw()

figure; SubjStatsSub(1).draw('tstat', [], 'q<0.05')
SubjStatsSub(1).variables



% load SubjStats if they already exist
load('Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Analysis_NIRS\gonogo\SubjStats-gonogo-2022-12-13.mat');

%% Save subject level stats to a folder
% Create one file with all subject level statistics for separate plots. 

% TT = [];
% for i = 1:1:numel(SubjStats)
%     T = SubjStats(i).table();
%     % add subject id 
%     subject = str2num(SubjStats(i).demographics.subject);
%     subject = repmat(subject, height(SubjStats(i).table()),1);
%     T = addvars(T,subject,'Before', "source");
% 
%     % add load 
%     load = {SubjStats(i).demographics.load};
%     load = repmat(load, height(SubjStats(i).table()),1);
%     T = addvars(T,load,'Before', "source");
% 
%     %add prepost
%     timing = {SubjStats(i).demographics.timing};
%     timing = repmat(timing, height(SubjStats(i).table()),1);
%     T = addvars(T,timing,'Before', "source");
% 
%     % append to large table
%     TT = vertcat(TT,T);
% end
% 
% 
% % save the stats tables
% stats_folder = ['C:\Users\henri\OneDrive\cibsr\BrexFit\figs\gonogo'];  % need to make this folder if it doesn't already exist
% filename = append('\2022-12-14-SubjStats.csv');
% writetable(TT, [stats_folder, filename]);

%% Assess subject level leverage for group model
% QA for Subject Level Stats
% generates a table with leverage for subj, condition & channel
groupLeverage = nirs.util.grouplevelleveragestats(SubjStats);

head(groupLeverage)

idx = find(groupLeverage.pval_subject < 0.05); %0.022
disp(groupLeverage(idx, :)); 

%TODO: Decide whether to exclude any subjects based on this 
% https://en.wikipedia.org/wiki/Leverage_(statistics) 
% 'Although an influential point will typically have high leverage, a high 
% leverage point is not necessarily an influential point.' 

% Remove outlier/highly influential subjects
j = nirs.modules.RemoveOutlierSubjects();
j.formula
j.cutoff = 0.05; 
GoodSubjs = j.run(SubjStats);
%% Change subject stats: Subtract go from nogo. 

% I have a question about how to combine two conditions in the subject
% level statistics to one condition [in nirs-toolbox]. Our experiment
% design consists of pre and post fNIRS measurements of a Go/Nogo task,
% with an intervention between. We are comparing two interventional
% conditions to a control condition (i.e., three levels). On the subject
% level we currently have beta estimates for both Go and Nogo trials. When
% I take all of these variables into a group level model I end up with a
% three way interaction effects model that’s very hard to interpret, and I
% don’t think it actually captures what we are trying to measure. See
% below: Cond = Go/Nogo, timing = pre/post, load =
% intervention1/intervention2/control-condition. j =
% nirs.modules.MixedEffects( ); j.formula = 'beta ~ -1 + cond*load*timing +
% (1|subject)';
% 
% I would like to convert the two beta estimates for Go and Nogo trials to
% one beta value that captures the Nogo-go contras, i.e., subtracting
% Go-trials from Nogo trials to get one beta capturing hemodynamic activity
% associated with response inhibition). I believe this would allow us to
% create a group model that considers pre-post difference in hemodynamic
% activity and its interaction with the various interventions (which is
% what we want to investigate). See below: j = nirs.modules.MixedEffects(
% ); j.formula = 'beta ~ -1 + load*timing + (1|subject)';
% 
% My question is thus the following: How do I subtract the Go from the Nogo
% trials (on a per channel basis) in the SubjStats data structure? I am
% able to subtract the betas from each other, but I don’t know how I should
% handle the remaining variables in the SubjStats data structure (t-stat,
% se, covb, etc.) to ensure that the mixed effects model gets all of the
% variables it requires.
% 





%% 2-nd level stats (group level) 

% Dummy coding scheme decides "contrasts":
% https://se.mathworks.com/help/stats/dummy-indicator-variables.html#mw_b475bf3b-7207-4a73-924e-e7c138b89769
% We want effects coding 

% Three and two-way interaction effects between condition, load and timing,
% and main effects of all. 
j = nirs.modules.MixedEffects( );
%j.formula = 'beta ~ -1 + cond*load*timing + (1|subject)'; %(1|age)
j.formula = 'beta ~ -1 + load*timing + (1|subject)'; %(1|age) after calculating one beta for nogo-go
j.dummyCoding = 'effects';
% (in final analysis) change defaults to: run robust stats, include diagnostics, and print
% what's happening while the model runs
  j.robust = true;
  j.include_diagnostics=true;
j.verbose = true;
GroupStats = j.run( SubjStats);

GroupStats.probe.defaultdrawfcn = '10-20'
figure; GroupStats.draw('tstat', [-5 5], 'q < 0.05') %

GroupStats.draw()
% save figures
figs_folder = ['Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Analysis_NIRS\gonogo\GroupStats\10-20']; 
GroupStats.printAll('tstat', [-5 5], 'q < 0.05', figs_folder, 'tif');

% save the stats tables
stats_folder = figs_folder; % need to make this folder if it doesn't already exist
writetable(GroupStats.table(), [stats_folder, '\2022-12-14-GroupStats.csv']);

head(GroupStats.table(),20)


FStats = GroupStats2.ftest(c);
FStats.table()
FStats.draw(10,'q<0.05') 

% Look at HRF
GroupStats2.HRF
figure; GroupStats3.HRF.draw()

%% 2nd attempt at 2nd level stats

% Main effect of timing. Controlling for subject and age (i.e subject and
% age as random variable)
j = nirs.modules.MixedEffects( );
j.formula = 'beta ~ -1 + timing + (1|subject) + (1|age)';
%j.dummyCoding = 'effects'; %'full'; 
% Question: which dummy coding?
j.verbose = true;
GroupStatsTiming = j.run(SubjStats);
GroupStatsTiming.draw('tstat', [-5 5], 'q < 0.05')

c =  {'pre-post'}; 
ContrastStats = GroupStatsTiming.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05')

% Main effect of gonogo (sanity check). Subject and age as random
% variables. 
j = nirs.modules.MixedEffects( );
j.formula = 'beta ~ -1 + cond + (1|subject) + (1|age)';
j.dummyCoding = 'full';
j.verbose = true;
GroupStatsTarget = j.run( SubjStats);
GroupStatsTarget.draw('tstat', [-5 5], 'q < 0.05')
c =  {'nonTarget-target'}; 
ContrastStats = GroupStatsTarget.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05')

% j = nirs.modules.MixedEffects( );
% j.formula = 'beta ~ -1 + cond + (1|subject) + (1|age)';
% j.dummyCoding = 'full';
% j.robust = true;
% j.verbose = true;
% GroupStatsTarget = j.run( SubjStats);
% GroupStatsTarget.draw('tstat', [-5 5], 'q < 0.05')


% Two-way interaction with lover order effects
j = nirs.modules.MixedEffects( );
j.formula = 'beta ~ -1 + timing*cond + (1|subject) + (1|age)';
j.dummyCoding = 'effects';
j.robust = true;
j.verbose = true;
GroupStats = j.run(SubjStats);
GroupStats.draw('tstat', [-5 5], 'q < 0.05')

% j = nirs.modules.MixedEffects( );
% j.formula = 'beta ~ -1 + load + (1|subject) + (1|age)'
% j.dummyCoding = 'full';
% j.verbose = true;
% GroupStats1 = j.run( SubjStats);

%% Contrasts
% Target - non-target
GroupStats.conditions 
%     {'in'                        }
%     {'out'                       }
%     {'pre'                       }
%     {'pre:load_in'               }
%     {'pre:load_out'              }
%     {'target'                    }
%     {'target:load_in'            }
%     {'target:load_out'           }
%     {'target:timing_pre'         }
%     {'target:timing_pre:load_in' }
%     {'target:timing_pre:load_out'}

c =  {'target'}; % Main effect of condition 
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05')

c =  {'out-in'}; % Main effect of condition 
ContrastStats = GroupStats.ttest(c);
ContrastStats.probe.defaultdrawfcn = '10-20'
ContrastStats.draw('tstat',[-5, 5],'q<0.05')

FStats = GroupStats.ftest(c);
head(FStats.table())
FStats.draw(10,'q<0.05') 

c =  {'pre'}; % Main effect of timing -- uninteresting in itself 
FStats = GroupStats3.ftest(c);
head(FStats.table())
FStats.draw(10,'q<0.05') 

% Interaction effect:
% Compares target to non-target and pre to post for in and out session
% compared to qr, (maybe do an ftest here? not sure if f or t is more appropriate for interaction effects)
c =  {'target:timing_pre:load_in', 'target:timing_pre:load_out'};
ContrastStats = GroupStats3.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05')

% Nope I don't think its appropriate with an F test, as it does not seem to
% provide directionality, which is something we do want for interaction
% effects. 
FStats = GroupStats3.ftest(c);
head(FStats.table())
FStats.draw(10, 'q<0.05') 

% Compares target to non-target and pre to post, and the outdoor to the
% indoor session
c =  {'target:timing_pre:load_out-target:timing_pre:load_in'};
ContrastStats = GroupStats3.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05')
% ok so there seems to be an interaction effect here, with greater increase
% in hbo in the outdoor session compared to indoor session 

GroupStats3.HRF

%% Attempt at subtracting Non-target-target


head(SubjStats(1).table(),20)

beta = -37.399

beta =

  -37.3990

se = 41.186

se =

   41.1860

beta/se

ans =

   -0.9081



% source detector type and ShortSeparation columns must be the same, while
% cond is supposed to be different

T = SubjStats(1).table();
% cond must be different 


find(T.source == 1 & T.detector == 1 & strcmp(T.type,{'hbo'}) & T.ShortSeperation == false)

find(strcmp(T.type,{'hbo'}))

find(strcmp(T.cond,{'target'}))

% perhaps split into two tables first
nonTargetTable = T(find(strcmp(T.cond,{'nonTarget'})), :)
targetTable = T(find(strcmp(T.cond,{'target'})), :)

head(nonTargetTable)
head(targetTable)
% for these tables row numbers match, we can therefore subtract 

newTbl = nonTargetTable(:,["source", "detector", "type", "ShortSeperation"]);

newTbl.beta = nonTargetTable.beta - targetTable.beta;
newTbl.cond = repmat({"nonTarget-target"},height(newTbl),1);
head(newTbl)

SubjStats(1).beta

SubjStats(1).variables

strcmp(SubjStats(1).conditions, 'nonTarget')










