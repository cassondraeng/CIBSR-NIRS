%% Script that analyzes BrexFit nback data 
% To run this script you need nirs-toolbox and QT-nirs 

% Last updated: 2022-12-19
% Created by: Henrikke Dybvik
% Using MATLAB R2021b


% Questions and remainding tasks: 
% TODO: run a mixed model with load and timing and see if that is more
% aligned with our desired output. 
% TODO: Find out which mixed effects model(s) we want to run. (might want
% to include age and baseline executive function in model)
% TODO: (Connected to the one above) Decide on contrasts. 


% TODO: Add option to plot all QT-nirs results in one figure for individual
% participants. 


%   NEED INPUT: 
% TODO: Inspect (Nback tasks are incomplete    '3010') 
% could be the first half, could also be the second half
% Need notes from QC excel sheet to figure this out


%% load nback data 
% Uncomment/comment dirlist and rootdir below according to os:

% [For Windows] get directory of files when using windows
rootdir = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS';
dirlist = dir(fullfile(rootdir, '**\*_nback\202*')); % get folders with nback data

% dirlist = dir(fullfile(rootdir, '**\*_nback\202*')); % [use this in the new year] get folders with nback data

% [For Mac]
%The rootdir needs to match how you login to the server, which is currently
%via Cisco Anyconnect with the following path entered into the "Go"
%connection smb://dibshome/Projects/NIRS_Projects/SpecializedFoundation/BREXFIT_S3
%rootdir = '/Volumes/BREXFIT_S3/Data/NIRS/';
%dirlist = dir(fullfile(rootdir, '*/*/*_nback/202*')); % get folders with nback data

% View folders and paths
%dirlist.name
%dirlist.folder

raw = []; % declare empty data set
for i=1:1:length(dirlist) % iterate through dirlist 
    % Select path depending on os (mac vs windows)
    path = append(dirlist(i).folder, '\', dirlist(i).name); % [For windows] find full path of individual folders
    %path = append(dirlist(i).folder, '/', dirlist(i).name); % [For mac] find full path of individual folders
    
    filedesc = dirlist(i).name;
    exp1 = '[a-z_]+_nback'; %matches one or more letters and underscores followed by _nback
    experiment = char(regexpi(dirlist(i).folder, exp1, 'match'));
    exp = '\d{4}'; %matches five consecutive digits: i.e. participant id
    participant_id = char(regexp(dirlist(i).folder, exp, 'match'));
    
    % filter out files that should be handled individually
    if strcmp(participant_id,'3008') & startsWith(experiment, 'qr_post')  
        % qr post nback: two files: 004=5.2mb, 005=8.1mb -> File 005 is correct full data file
        if endsWith(filedesc, '005') % this is the correct file
            % load
            raw_ind = nirs.io.loadDirectory(path); % load individual folders/files
            raw = [raw raw_ind]; % append individual files to data set 
            %filedesc % correct file, remove the other one 
        else 
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end
    elseif strcmp(participant_id,'3016') & startsWith(experiment, 'qr_post')
        % two files: 005=227kb, 006=8.5, File 006 is correct full data file
        if endsWith(filedesc, '006')
            raw_ind = nirs.io.loadDirectory(path); % load individual folders/files
            raw = [raw raw_ind]; % append individual files to data set 
        else
            fprintf('NOT LOADING THIS FILE: %s+n', filedesc)
            continue
        end
    elseif strcmp(participant_id,'3014') & startsWith(experiment, 'out_post')
        % two files: 003=839kb, 004=8mb, File 004 is correct full data file
        if endsWith(filedesc, '004')
            raw_ind = nirs.io.loadDirectory(path); % load individual folders/files
            raw = [raw raw_ind]; % append individual files to data set 
        else
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end
    else
        try
            raw_ind = nirs.io.loadDirectory(path); % load individual folders/files
            raw = [raw raw_ind]; % append individual files to data set 
        catch expression %expression is an MException struct
            disp(participant_id); disp('caused error')
            fprintf(1,'The identifier was:\n%s',expression.identifier);
            fprintf(1,'There was an error! The message was:\n%s',expression.message);
            continue
        end
   end
end


%% Visual inspection 
% (view entire dataset)
raw % view data set structure
raw.draw() % graph entire data set
figure; raw(2).draw() % graph single file

% There are files without stimui (i.e., triggers) here

demoTable = nirs.createDemographicsTable(raw); % create demographics table
disp(demoTable); % view demoghraphics table

% plot montage
figure; raw(1).probe.draw()

raw(1).probe.defaultdrawfcn='10-20';
figure; raw(1).probe.draw()

raw(1).probe.defaultdrawfcn='?;' % view possible draw configurations

raw(1).probe.defaultdrawfcn='3D mesh';
raw(1).probe.defaultdrawfcn='3D mesh (top)' ;


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
figure; raw(2).draw()

% Event markers: 
% one_back = 1
% two_back = 2
% zero_back = 4
% non_target = 8
% target = 16


% rename stimuli 
j = nirs.modules.RenameStims()
j.listOfChanges = {
    'stim_channel1', 'one_back';
    'stim_channel2', 'two_back';
    'stim_channel4', 'zero_back';
    'stim_channel8', 'non_target';
    'stim_channel16', 'target'
    };
raw = j.run(raw);

% 1, 2 and 4 are the start of a block, whereas 8 and 16 indicate target
% /non-target. 

% discard stimuli not in use for block analysis (for now)
j=nirs.modules.DiscardStims;
j.listOfStims={
     'non_target';
     'target'};
raw=j.run(raw);


% % for ERP analysis (to be removed from this script)
% job=nirs.modules.DiscardStims;
% job.listOfStims={
%      'one_back';
%      'two_back';
%      'zero_back'};
% rawerp=job.run(raw);
% figure;rawerp(2).draw()

% remove files without simuli 
j = nirs.modules.RemoveStimless();
raw = j.run(raw);

% Calculate block durations:
% loaded trial info from brexfit study 3 cognitive task script and copied
% the output her so that that the trial info file doesen't need to be
% loaded evey time 
% dur = trialInfo.trialNum * 2.35
dur =[
   18.8000
   21.1500
   23.5000
   25.8500
   28.2000
   25.8500
   25.8500
   21.1500
   25.8500
   25.8500
   28.2000
   23.5000
   28.2000
   23.5000
   23.5000]; % duration vector

N = [ 0
     2
     1
     0
     1
     2
     0
     2
     1
     0
     1
     2
     0
     2
     1]; %trialInfo.N


% find index and subsquently durations for zero, one and two back blocks:
idx_0back = find(N==0); % find index of zero-back blocks
dur_0back = dur(idx_0back); % select parts of duration vector that corresponds to zero-back

idx_1back = find(N==1);
dur_1back = dur(idx_1back);

idx_2back = find(N==2);
dur_2back = dur(idx_2back);


% change stimuli duration (must rename and remove stimless files first to make it work)
for i=1:1:numel(raw) % set duration for 0-back
    try
        raw(i).stimulus.zero_back.dur(:)=dur_0back;
    catch expression %expression is an MException struct
        fprintf(1,'raw(%i) caused and error when setting duration for 0-back \n',i)
        fprintf(1,'The identifier was:\n%s\n',expression.identifier);
        fprintf(1,'The message was:\n%s\n',expression.message);
        continue
    end
end

for i=1:1:numel(raw) % set duration for 1-back
    try
        raw(i).stimulus.one_back.dur(:)=dur_1back;
    catch expression %expression is an MException struct
        fprintf(1,'raw(%i) caused and error when setting duration for 1-back  \n',i)
        fprintf(1,'The identifier was:\n%s\n',expression.identifier);
        fprintf(1,'The message was:\n%s\n',expression.message);
        continue
    end
end

for i=1:1:numel(raw) % set duration for 2-back
    try
        raw(i).stimulus.two_back.dur(:)=dur_2back;
    catch expression %expression is an MException struct
        fprintf(1,'raw(%i) caused an error when setting duration for 2-back  \n',i)
        fprintf(1,'The identifier was:\n%s\n',expression.identifier);
        fprintf(1,'The message was:\n%s\n',expression.message);
        continue
    end
end
% Look at the output and handle any errors individually

% handle errors individually 
figure; raw(3).draw()
raw(3).stimulus.zero_back.dur(:)=dur_0back(2:end);
figure; raw(3).draw()

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% TODO Figure out what's going on here and correct stimuli duration.
j = 32;
figure; raw(j).draw()
raw(j).demographics.subject % what is going on here ? The file is incomplete
raw(j).demographics.experiment
raw(j).stimulus.one_back.count
raw(j).stimulus.two_back.count
raw(j).stimulus.zero_back.count
% the sequence for i = 29 is 
% 0 
% 2
% 1
% 0
% 1
% 2
% 0
% 2
N = [ 0
     2
     1
     0
     1
     2
     0
     2
     1
     0
     1
     2
     0
     2
     1]; 
% 0 
% 2
% 1
% 0
% 1
% 2
% 0
% 2
% could be the first half, could also be the second half
% Need notes from QC excel sheet to figure this out

% visual inspection of individual 
figure; raw(4).draw()

figure; raw.draw()

%% rename variables according to condition. (pre/post) + (qr/in/out)
% We have the variable timing (pre, post) and the variable load (qr/out/in)
% These can be placed in demographics table and be used as predictors or 
% covariates in the group level model later. 
% convert subject id to numeric to be able to add demmogrpahic info from 
% .csv file later, and create a character version of subject_id in case
% that's neccessary for mixed model later. 

for i =1:numel(raw)
    % rename experiment variable 
    newExpStr = erase(raw(i).demographics.experiment, '_nback');
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
    % subject id
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

% Maybe save excel sheet as .csv
% add path for mac
additional_demographics = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\Demographics\BrainResearchExercis-Demographics2.xlsx';
j.demoTable = readtable(additional_demographics);

% We are going to match the subject column in the above table
j.varToMatch = 'subject'; 
% Note "subject" (or whatever you are matching based on) needs to be an
% entry in BOTH the CSV file and the nirs.core.<> class that you are
% loading the demographics to. "Subject" need to be numerical values for
% this to work.

raw = j.run(raw); % Can run it on SubjStats as well. 

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


%% Order dataset so post is listed first (and will be used as reference group)

for i=1:numel(raw)
    if strcmp(raw(i).demographics.timing,'post')
        % identify first occurence of post
        rawpostfirst = raw(i);
        idx_remove = i;
        break
    else
        continue
    end
end
raw(idx_remove) = [];
raw = [rawpostfirst raw];

%% label the short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15;
raw = j.run(raw);

% check that short-channels are labeld (command might produce error if the
% field ShortSeparation does not exist
raw(1).probe.link.ShortSeperation == 1 
%% Quality check 

j = nirs.modules.QT()
j.qThreshold =0.7 % what is q??
j.sciThreshold = 0.6
j.fCut = [0.5 2.0]
ScansQuality = j.run(raw);
% ScansQuality = j.run(raw(1:5)) % if we want to run on a subset

% view group results
ScansQuality.drawGroup('sq');
ScansQuality.drawGroup('sci');
ScansQuality.drawGroup('psp');
ScansQuality.drawGroup('bar');

% view individual results
for ii = 1:2 %numel(ScansQuality)
   ScansQuality(ii).drawGroup('sq');
   ScansQuality(ii).drawGroup('sci');
   ScansQuality(ii).drawGroup('psp');
   disp(ii)
   disp(raw(ii).demographics)
end

%% Identify how many channels are bad (percentage) and where they are
% located
% Its possible to calculate survivalrate at a couple thresholds 

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

%% decide how much data to kick out at this stage
% look at both sci and psp
% regular cuttoff at psp is 0.1
% for sci I'm tinking 0.5 (or lower, depending on how much data we're 
% loosing ) or 0.6 at max



%% preprocessing for block analysis

% basic preproccesing pipeline 
% j = nirs.modules.default_modules.basic_preprocessing();
% list = nirs.modules.pipelineToList(j);
% disp(list);
% %  {1×1 nirs.modules.ImportData    }
% %     {1×1 nirs.modules.RemoveStimless}
% %     {1×1 nirs.modules.FixNaNs       }
% %     {1×1 nirs.modules.Resample      }
% %     {1×1 nirs.modules.OpticalDensity}
% %     {1×1 nirs.modules.BeerLambertLaw}
% %     {1×1 nirs.modules.TrimBaseline  }
% %     {1×1 nirs.modules.ExportData    }
% list{4}.Fs = 1; %downsample to 1 Hz for now to speed things up
% list{7}; % Trim baseline, leave defaults of 30 s pre/post stimuli
% j = nirs.modules.listToPipeline(list) % convert back to pipeline
% hb = j.run(raw);


% Alternative
j = nirs.modules.FixNaNs();
j = nirs.modules.Resample(j);
j.Fs = 5;
j = nirs.modules.OpticalDensity(j);
j = nirs.modules.BeerLambertLaw(j);
hb = j.run(raw);

% visual inspection 
figure; hb(1).draw

% Label short separation channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15;
hb = j.run(hb);

% check that short-channels are labeld (command might produce error if the
% field ShortSeparation does not exist
hb(1).probe.link.ShortSeperation == 1 
%% 1st level stats (Subject Level) w/Short Channel regression 

load('Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Analysis_NIRS\preprocessed\Hb-2022-12-02.mat');
% Ran this on Hb-2022-12-02.mat. Startet 09:19. Fs = 5Hz. Let's see how
% long it takes. one file took about 5 min (completed 09:24). 4 files took
% 11 min (completed 09:31). Done by 12:00, bu dont knwo when
j = nirs.modules.GLM();  
disp(j) % show default parameters of GLM, use AR-IRLS
j.verbose = true; % add verbose function (to show progress)
j.AddShortSepRegressors = 1; %add short channel regressors
SubjStatsSS = j.run(hb); % run job and save the output to SubjStatsSS variable


load('Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Analysis_NIRS\preprocessed\SubjStatsSS-2022-12-03.mat');

head(SubjStatsSS(1).table(),20)
%% leverage at subject level
% QA for Subject Level Stats

% generates a table with leverage for subj, condition & channel
groupLeverage = nirs.util.grouplevelleveragestats(SubjStatsSS);

%head(groupLeverage)
% find high leverage subjects
idx = find(groupLeverage.pval_subject < 0.05); 
disp(groupLeverage(idx, :)); 

% find high leverage conditions
idx = find(groupLeverage.pval_condition < 0.05); 
disp(groupLeverage(idx, :)); 


%TODO: Decide whether to exclude any subjects based on this 
% https://en.wikipedia.org/wiki/Leverage_(statistics) 
% 'Although an influential point will typically have high leverage, a high 
% leverage point is not necessarily an influential point.' 

% Remove outlier/highly influential subjects
j = nirs.modules.RemoveOutlierSubjects();
j.formula
j.cutoff = 0.05; % 0.051; % Huppert: p-value cutoff (2 subjects almost exactly at p = 0.05, in original analysis
% with v 615 were removed at 0.05, current version set to 0.051 to exclude them)
GoodSubjs = j.run(SubjStatsSS);

%% 2nd (group) level stats
% Maybe run several models ?
% We have three predictors do we not? 

% We are interested in the following effects:
% - main effect of condition (i.e. high and low load  (qr, in, out))
% - main effect of timing (pre-post-change) (regardless of condition, is there an effect of anyhting???)
% - interaction effect, is there a differential effect from pre-post based
%   on condition 


% Dummy coding scheme decides "contrasts":
% https://se.mathworks.com/help/stats/dummy-indicator-variables.html#mw_b475bf3b-7207-4a73-924e-e7c138b89769
% We want effects coding 

% SubjStatsSS.conditions
% {'one_back' }
% {'two_back' }
% {'zero_back'} % last category so this must be the one that's subtracted

% This must apply for load and timing as well. 

% SubjStatsSS(1).demographics.load
% ans =
%     'in'
% SubjStatsSS(1).demographics.timing
% ans =
%     'pre'
% SubjStatsSS(2).demographics.load
% ans =
%     'out'



% 1st model: (the one we want)
load('Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Analysis_NIRS\preprocessed\GroupStats-2022-12-06.mat');
% Three-way interaction with all lower order effects:
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + load*timing*cond + (1|subject)'; % include covariates here (might want
% to include age and baseline executive function in model)
j.dummyCoding = 'effects'; % Cannot use 'full' here
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStats = j.run(SubjStatsSS);
% Started 12:38. Let's see how long it takes. First iteration finished
% 13:00. finished a bit past 1400, like 1408. But was not saved. 

% Conditions relative to baseline (Beer-Lambert law baseline)
% Draw the probe using t-stat & false discovery rate (q value)
GroupStats.draw('tstat', [], 'q < 0.05')
head(GroupStats.table(),20)



% 2nd model (not the one we want)
load('Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Analysis_NIRS\preprocessed\GroupStats-2022-12-03.mat');
% Three-way interaction with all lower order effects: (we're not estimating
% the intercept)
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + load*timing*cond + (1|subject)'; % include covariates here
j.dummyCoding = 'reference';
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStats2 = j.run(SubjStatsSS);
% Finished solving: time elapsed 1667.4378s ( = 27 min)
GroupStats2.draw('tstat', [], 'q < 0.05')
head(GroupStats2.table(),20)

% source    detector     type                     cond                     beta        se        tstat      dfe        p           q        minDiscoverableChange    RelativePower
%     ______    ________    _______    __________________________________    ________    _______    ________    ___    _________    ________    _____________________    _____________
% 
%       1          1        {'hbo'}    {'two_back'                      }      3.6509     1.3853      2.6354    211    0.0090278      0.0734           3.4571               0.31002   
%       1          1        {'hbo'}    {'zero_back'                     }     -2.5012     1.3697     -1.8261    211     0.069249     0.24645           3.4182               0.31355   
%       1          1        {'hbo'}    {'post'                          }      2.6518     1.9198      1.3813    211      0.16865     0.41025           4.7908               0.22371   
%       1          1        {'hbo'}    {'out'                           }     -4.6764     1.6994     -2.7518    211    0.0064429    0.059132           4.2409               0.25272   
%       1          1        {'hbo'}    {'qr'                            }     -3.3009     1.3955     -2.3653    211      0.01892     0.11424           3.4826               0.30776   
%       1          1        {'hbo'}    {'two_back:timing_post'          }     -5.4412     2.9601     -1.8382    211     0.067441     0.24294           7.3869               0.14509   
%       1          1        {'hbo'}    {'zero_back:timing_post'         }    -0.50448     2.9591    -0.17049    211      0.86479     0.95021           7.3844               0.14514   
%       1          1        {'hbo'}    {'two_back:load_out'             }    -0.28752     2.6702    -0.10768    211      0.91436     0.96738           6.6636               0.16084   
%       1          1        {'hbo'}    {'zero_back:load_out'            }      6.3774     2.6672       2.391    211     0.017679     0.10876            6.656               0.16102   
%       1          1        {'hbo'}    {'two_back:load_qr'              }     -3.5453     2.3081     -1.5361    211      0.12602     0.34926           5.7598               0.18608   
%       1          1        {'hbo'}    {'zero_back:load_qr'             }      2.5939     2.2878      1.1338    211      0.25815     0.51835           5.7091               0.18773   
%       1          1        {'hbo'}    {'post:load_out'                 }      3.2226     2.9929      1.0768    211      0.28282     0.54747           7.4688                0.1435   
%       1          1        {'hbo'}    {'post:load_qr'                  }     -3.1665     2.8514     -1.1105    211      0.26806     0.52875           7.1158               0.15062   
%       1          1        {'hbo'}    {'two_back:timing_post:load_out' }      1.7313     4.3446     0.39848    211      0.69068     0.86682           10.842              0.098853   
%       1          1        {'hbo'}    {'zero_back:timing_post:load_out'}     -2.6207     4.3426    -0.60349    211      0.54683     0.78206           10.837              0.098899   
%       1          1        {'hbo'}    {'two_back:timing_post:load_qr'  }      2.7644     4.1218     0.67068    211      0.50316     0.75215           10.286                0.1042   
%       1          1        {'hbo'}    {'zero_back:timing_post:load_qr' }      3.6781     4.1233     0.89205    211      0.37338     0.63661            10.29               0.10416   
%       1          1        {'hbr'}    {'two_back'                      }      2.1423    0.66347       3.229    211    0.0014408    0.023176           1.6557               0.64733   
%       1          1        {'hbr'}    {'zero_back'                     }      1.4283    0.65629      2.1764    211     0.030638     0.15511           1.6378               0.65441   
%       1          1        {'hbr'}    {'post'                          }     0.48693     0.7241     0.67245    211      0.50203     0.75215            1.807               0.59312   


% (in final analysis) change defaults to: run robust stats, include diagnostics, and print
% what's happening while the model runs
%  j.robust = true;
%  j.include_diagnostics=true;
%  j.verbose = true;

% test on restructured dataset 

j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + load*timing*cond + (1|subject)';
j.dummyCoding = 'effects'; % Cannot use 'full' here
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStatsRes = j.run(SubjStatsSS2);
%% 2nd attempt at group level stats (pre-post) 

% Two-way interaction with all lower order effects:
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + timing*cond + (1|subject) + (1|age)'; % might want to include baseline executive function in model)
j.dummyCoding = 'effects';  % cannot use full here
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStats = j.run(SubjStatsSS);

GroupStats.draw('tstat', [], 'q < 0.05')
% There is one t value thats 10^12


% Two-way interaction with all lower order effects:
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + timing*cond + (1|subject) + (1|age)'; % might want to include baseline executive function in model)
j.dummyCoding = 'effects';  % cannot use full here
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStats2 = j.run(GoodSubjs);
GroupStats2.draw('tstat', [-10 10], 'q < 0.05')
head(GroupStats2.table(),20)


% Main effects of timing 
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + timing + (1|subject) + (1|age)'; % might want to include baseline executive function in model)
j.dummyCoding = 'full';  
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStats3 = j.run(GoodSubjs);
head(GroupStats3.table(),20)
GroupStats3.draw('tstat', [-10 10], 'q < 0.05')

% Main effects of nback 
j = nirs.modules.MixedEffects();
j.formula = 'beta ~ -1 + cond + (1|subject) + (1|age)'; % might want to include baseline executive function in model)
j.dummyCoding = 'full';  
%j.robust = true;
%j.include_diagnostics=true;
j.verbose = true;
GroupStats3 = j.run(GoodSubjs);
head(GroupStats3.table(),20)
GroupStats3.draw('tstat', [-10 10], 'q < 0.05')

%% contrasts
% 1st model
% disp(GroupStats.conditions);
%   {'in'                          }
%     {'one_back'                    }
%     {'one_back:load_in'            }
%     {'one_back:load_out'           }
%     {'one_back:timing_pre'         }
%     {'one_back:timing_pre:load_in' }
%     {'one_back:timing_pre:load_out'}
%     {'out'                         }
%     {'pre'                         }
%     {'pre:load_in'                 }
%     {'pre:load_out'                }
%     {'two_back'                    }
%     {'two_back:load_in'            }
%     {'two_back:load_out'           }
%     {'two_back:timing_pre'         }
%     {'two_back:timing_pre:load_in' }
%     {'two_back:timing_pre:load_out'}

GroupStats.conditions

% -------------------Main effects:--------------------------------------
% Define contrast(s) for n-back:
% 2back-0back, 2back-1back, 1back-0back
c = {'two_back', 'two_back-one_back', 'one_back'};
ContrastStats = GroupStats.ttest(c); % Calculate stats with the ttest function
ContrastStats.draw('tstat',[-5, 5],'q<0.05')

% Define contrasts for load (only): 
% Fixed effect (Main effect of load)
% out-qr, in-qr, out-in
c =  {'out', 'in','out-in'};
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05')
% FStats = GroupStats.ftest(c); % p-values are the same, depends on if you
% want an F-statistic or t-statistic. 
% FStats.draw(10,'q<0.05') 

% Define contrasts for timing
% Main Effect of timing
% post - pre
c =  {'-pre'};
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05') 

% --------------Two-way Interaction effects-------------------------------

%     {'one_back:load_in'            }
%     {'one_back:load_out'           }

%     {'two_back:load_in'            }
%     {'two_back:load_out'           }

%     {'one_back:timing_pre'         }
%     {'two_back:timing_pre'         }

%     {'pre:load_in'                 }
%     {'pre:load_out'                }



% pre-post difference interaction effect with load
c= {'pre:load_in', 'pre:load_out'}
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05') 


% pre-post difference interaction effect with load
c= {'pre:load_out-pre:load_in'}
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05') 

% --------------Three-Way Interaction Effects-----------------------------

%     {'one_back:timing_pre:load_in' }
%     {'one_back:timing_pre:load_out'}

%     {'two_back:timing_pre:load_in' }
%     {'two_back:timing_pre:load_out'}


% (1back-0back)pre-post(in-qr), (1back-0back)pre-post(out-qr)
c= {'one_back:timing_pre:load_in', 'one_back:timing_pre:load_out'};
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05') 


% (2back-0back)pre-post(in-qr), (2back-0back)pre-post(out-qr)
c= {'two_back:timing_pre:load_in','two_back:timing_pre:load_out'}; % not sure we can simply reverse signs here. 
ContrastStats = GroupStats.ttest(c);
ContrastStats.draw('tstat',[-5, 5],'q<0.05') 

head(ContrastStats.table)
idx = find(ContrastStats.q < 0.05); 
SignificantContrastStats = ContrastStats.table; 
SignificantContrastStats = SignificantContrastStats(idx, :); 
disp(SignificantContrastStats);

%ContrastStats.draw('tstat',[-10 10],'q<0.05')
ContrastStats.draw('tstat',[],'q<0.05')

% =================================
% model

disp(GroupStats2.conditions);
% define contrast(s) 
c = {'two_back-one_back'}

% ===================================
%  model
disp(GroupStats3.conditions);


