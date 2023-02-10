%% Script that analyzes BrexFit during data
% and combines nirs data files, e.g., outdoor during files

% To run this script you need nirs-toolbox and QT-nirs 
% Last updated: 2022-12-12
% Created by: Henrikke Dybvik
% Using MATLAB R2021b



% TODO

% combine during files for qr, in, out into one dataset
% quality check
% Decide how to handle shorter recordings (less than 20 min) and do that. 
% Compare various preprocessing pipelines for fc 

% visual insepction of 
        % 3004 NIRS device disconnected during warmup -> reconnect at 12:30,
        % Trigger placed at end of exercise & beginning, The file with the trigger
        % is 003 and has the correct desc.json file, and  there is a mis_nirs_file
        % labeled 004…but neither appear to have the right .nirs file size

        % 3016 two files: 003=2.2mb, 004=19mb *See note: 1) Bike trainer was not
        % plugged in so garmin paired to power meter 2) Had to stop files/data
        % collection at 2:21 into the warm-up 3) File #2 on phone & garmin have the
        % exercise data from the start of 2:21 onward 4) NIRS files also has 2
        % files; 1 for the start of exercise to 2:21 and 1 file for exercise from
        % 2:21 onward…tiggers are included 3 times (1 for restart, 1 for end of
        % warmup/start of target HR zone, and one for end of exercise

        % HD: Perhaps it's best to load both files and combine them, then
        % decide which time period to use? 





% Questions:
% - I assume we should use the same preprocessing pipeline for qr, in and
% out? (and not individual preprocessing pipelines) Have written the code
% as if we're using one preproccesing pipeline. 

% - I think the group level connectivity breaks in part due to different
% montages being used. 

% https://www.sciencedirect.com/science/article/pii/S0304394019307074
% First, we used a threshold for signal to noise ratio (SNR) of 6.67 (∼ 15%
% coefficent of variation (CV); SNR = 1/CV * 100, HOMER2 PruneChannels
% function). With this procedure 10.0% of the channels were regarded as too
% noisy and not included in the further analysis


%% Notes
% Has been done:
% Combine during files for out

% https://github.com/Artinis-Medical-Systems-B-V/SignalQualityIndex
% https://support.nirx.de/archives/knowledge/can-i-merge-two-files-into-one
%% 1. loading QR during

rootdir = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS';
qr_dirlist = dir(fullfile(rootdir, '**\qr_during\2022*')); % get folders (this search is what takes time)

%qr_dirlist.name % gives the name of the folder
%qr_dirlist.folder % gives the path to the folder
% --> must combine the two above to get the folder path

raw_qr = []; % declare empty data set
for i=1:1:length(qr_dirlist) % iterate through dirlist 
    path = append(qr_dirlist(i).folder, '\', qr_dirlist(i).name); % find full path of individual folders

    filedesc = qr_dirlist(i).name;
    exp1 = '[a-z_]+_during'; %matches one or more letters and underscores followed by _nback
    experiment = char(regexpi(qr_dirlist(i).folder, exp1, 'match'));
    exp = '\d{4}'; %matches five consecutive digits: i.e. participant id
    participant_id = char(regexp(qr_dirlist(i).folder, exp, 'match'));

    if strcmp(participant_id,'3016') & startsWith(experiment, 'qr_during')  
        % 3016 qr_during two files: 003=2.5mb, 004=15.5mb, File 004 is correct full data file
        if endsWith(filedesc, '004') % this is the correct file
            % load
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_qr = [raw_qr raw]; % append individual files to data set 
            %filedesc % correct file, remove the other one 
        else 
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end    
    elseif strcmp(participant_id,'3023') & startsWith(experiment, 'qr_during')  
        % 3023 3 files due to restarts, File 005 has the most data and is from time
        % 15:30 until the end of the video: Device stopped working a few minutes 
        % into the video, but we kept the video running while we tried to restart. 
        % There are multiple nirs files and the “restart” (i.e., file 005 is from 
        % 15:30 until the end of the video)
        if endsWith(filedesc, '005') % this is the correct file
            % load
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_qr = [raw_qr raw]; % append individual files to data set 
            %filedesc % correct file, remove the other one 
        else 
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end
    else
        try
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_qr = [raw_qr raw]; % append individual files to data set 
        catch expression %expression is an MException struct
            disp(participant_id); disp('caused error')
            fprintf(1,'The identifier was:\n%s',expression.identifier);
            fprintf(1,'There was an error! The message was:\n%s',expression.message);
            continue
        end
    end    
end

% 3030 has two files, I guess qr_during2 is correct?

%raw_qr.draw() % visual inspection 

demoTable_qr = nirs.createDemographicsTable(raw_qr); % view 
disp(demoTable_qr);

raw_qr(1).demographics.experiment


%% correct stimuli name and duration [qr during]

raw_qr.stimulus

% discard aux stimuli
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_aux1'; 
    'stim_aux2';
    'stim_aux3';
    'stim_aux4';
    };
raw_qr=j.run(raw_qr);

% we have stim channel 1 and stim channel 2

% Add stimuli duration based on last timestamp: 
% (MIGHT NOT BE THE APPROACH WE WANT TO USE)
% last timestamp minus 5 minutes and then fifteen minutes before that. So
% 20 minutes from the last timestamp, and 15 min duration. Except if the
% duration is shorter than expected (25 min), due to it ending prematurely.
% create a new stimuli marker if there isn't one

% it might be easier to erase all markers and just add 

for i= 1:1:length(raw_qr)
    last_timestamp = raw_qr(i).time(end);
    if last_timestamp > 20*60
        % subtract 20 minutes from last timestamp and set that as onset
        onset = last_timestamp - 20*60;
        dur = 60*15;
        name = 'qr_during';
        amp = 1;
        qr_during = nirs.design.StimulusEvents(name, onset, dur, amp);
        % add stimuli event to nirs file
        raw_qr(i).stimulus('qr_during') = qr_during;
        disp('recording longer than 20 min ')
    else
        % if the recording is less than 20 minutes long, 
        % skipp it for now. I assume these must be handled individually
        % based on QC notes. 
        disp('skipping short session')
        continue
        % TODO: Decide how long these stimuli markers should be
    end
end

figure; raw_qr(2).draw()

% discard stimuli not in use
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_channel1'; 
    };
raw_qr=j.run(raw_qr);


% label the short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15; % check units 
raw_qr = j.run(raw_qr);

% check that short-channels are labeld (command might produce error if the
% field ShortSeparation does not exist
raw_qr(3).probe.link.ShortSeperation == 1 % some 1s should appear

%% Signal quality check via QT-nirs [QR]

% Extract task period and evaluate quality metrics over task period only
raw_qr_ext = extract_task_period(raw_qr);

j = nirs.modules.QT();
j.qThreshold = 0.6;
j.sciThreshold = 0.5;  % This paper: 10.1117/1.NPh.9.1.015001 used 0.5 as threshold
%j.pspThreshold = 0; %.5; % check with psp threshold = 0 
j.fCut = [0.5 3]; % corresponds to HR between 30 and 180
ScansQuality_qr = j.run(raw_qr_ext);

%ScansQuality.probe.defaultdrawfcn = '10-20';

% View group results
ScansQuality_qr.drawGroup('sq');
ScansQuality_qr.drawGroup('sci');
ScansQuality_qr.drawGroup('psp');
ScansQuality_qr.drawGroup('bar');
ScansQuality_qr.drawGroup('sqmask');

% View individual results 
for ii = 1:numel(ScansQuality_qr)
   %figure(ii),
   ScansQuality_qr(ii).drawGroup('sq');
   ScansQuality_qr(ii).drawGroup('sci'); 
   ScansQuality_qr(ii).drawGroup('psp');
end

%ScansQuality.qMats
%ScansQuality.probe

% Identify how many channels are bad (percentage) and where they are
% located

% Declare variables that will hold overall survival rate
ChTotALL = 0;
badChTotALL = 0;
ChLongALL = 0;
badChLongALL = 0;
ChssALL = 0;
badChssALL = 0;


fprintf('Survival rate:\n');
for i=1:1:numel(ScansQuality_qr)
    % find total number of channels
    Ch = raw_qr_ext(i).probe.link;
    ChTot = height(Ch)/2; % divide by two to account for hbo/hbr
    
    % find bad channels
    idxBadCh = find(ScansQuality_qr(i).qMats.MeasListAct==0);
    badChTot = length(idxBadCh)/2;
    
    % find long channels
    idxChLong = find(raw_qr_ext(i).probe.link.ShortSeperation == 0);
    ChLong = length(idxChLong)/2;
    
    % find bad long channels
    idxBadChLong = find((ScansQuality_qr(i).qMats.MeasListAct==0) & (raw_qr_ext(i).probe.link.ShortSeperation == 0));
    badChLong = length(idxBadChLong)/2;
    
    % find short channels
    idxChss = find(raw_qr_ext(i).probe.link.ShortSeperation == 1);
    Chss = length(idxChss)/2;
    
    % find bad short channels
    idxBadChss = find((ScansQuality_qr(i).qMats.MeasListAct==0) & (raw_qr_ext(i).probe.link.ShortSeperation == 1));
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

%% preprocessing qr



%% 2. loading IN during
% Loading indoor during data using the same procedure as above (1.) 

rootdir = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS';
in_dirlist = dir(fullfile(rootdir, '**\in_during\2022*')); % get folders

raw_in = []; % declare empty data set
for i=1:1:length(in_dirlist) % iterate through dirlist 
    path = append(in_dirlist(i).folder, '\', in_dirlist(i).name); % find full path of individual files
    
    filedesc = in_dirlist(i).name;
    exp1 = '[a-z_]+_during'; %matches one or more letters and underscores followed by _nback
    experiment = char(regexpi(in_dirlist(i).folder, exp1, 'match'));
    exp = '\d{4}'; %matches five consecutive digits: i.e. participant id
    participant_id = char(regexp(in_dirlist(i).folder, exp, 'match'));

    if strcmp(participant_id,'3004') & startsWith(experiment, 'in_during')  
        % 3004 NIRS device disconnected during warmup -> reconnect at 12:30,
        % Trigger placed at end of exercise & beginning, The file with the trigger
        % is 003 and has the correct desc.json file, and  there is a mis_nirs_file
        % labeled 004…but neither appear to have the right .nirs file size

        if endsWith(filedesc, '003') % load 003 file for now. Can be changed. Might want to visually inspect both. 
            % load
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_in = [raw_in raw]; % append individual files to data set 
        else 
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end    
    elseif strcmp(participant_id,'3007') & startsWith(experiment, 'in_during')  
        % 3007 two files: 003=2mb, 004=19.7 - No notes, but file 004 has trigger
        % file and looks the right size…also have misc_nirs_file but am unsure what
        % happened to create that
        if endsWith(filedesc, '004') %load 004
            % load
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_in = [raw_in raw]; % append individual files to data set 
        else 
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end
    elseif strcmp(participant_id,'3016') & startsWith(experiment, 'in_during')  
        % 3016 two files: 003=2.2mb, 004=19mb *See note: 1) Bike trainer was not
        % plugged in so garmin paired to power meter 2) Had to stop files/data
        % collection at 2:21 into the warm-up 3) File #2 on phone & garmin have the
        % exercise data from the start of 2:21 onward 4) NIRS files also has 2
        % files; 1 for the start of exercise to 2:21 and 1 file for exercise from
        % 2:21 onward…tiggers are included 3 times (1 for restart, 1 for end of
        % warmup/start of target HR zone, and one for end of exercise

        % HD: Perhaps it's best to load both files and combine them, then
        % decide which time period to use? 
        if endsWith(filedesc, '004') % load for 
            % load
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_in = [raw_in raw]; % append individual files to data set 
        else 
            fprintf('NOT LOADING THIS FILE: %s\n', filedesc)
            continue
        end
    else
        try
            raw = nirs.io.loadDirectory(path); % load individual folders/files
            raw_in = [raw_in raw]; % append individual files to data set 
        catch expression %expression is an MException struct
            disp(participant_id); disp('caused error')
            fprintf(1,'The identifier was:\n%s',expression.identifier);
            fprintf(1,'There was an error! The message was:\n%s',expression.message);
            continue
        end
    end    
end

% There was an error reading one file (participant 3004). 

% check if all data was loaded and flag if not
if isequal(length(in_dirlist), length(raw_in))
    disp('files might all be ok');
else
    disp('not all files seen to have been loaded correctly');
end



demoTable_in = nirs.createDemographicsTable(raw_in);
disp(demoTable_in);


raw_in(10).demographics.subject
raw_in(10).demographics.experiment

raw_in(11).demographics.subject
raw_in(11).demographics.experiment
raw_in(11).demographics.remarks

%% correct stimuli name and duration [IN during]

% discard aux stimuli
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_aux1'; 
    'stim_aux2';
    'stim_aux3';
    'stim_aux4';
    };
raw_in=j.run(raw_in);


% Add stimuli markers based on the end of the file
for i= 1:1:length(raw_in)
    last_timestamp = raw_in(i).time(end);
    if last_timestamp > 20*60
        % subtract 20 minutes from last timestamp and set that as onset
        onset = last_timestamp - 20*60;
        dur = 60*15;
        name = 'out_during';
        amp = 1;
        in_during = nirs.design.StimulusEvents(name, onset, dur, amp);
        % add stimuli event to nirs file
        raw_in(i).stimulus('in_during') = in_during;
        disp('recording longer than 20 min')
    else
        % if the recording is less than 20 minutes long, 
        % skipp it 
        disp('skipping short session')
        continue
        % TODO: Decide how long these stimuli markers should be
    end
end

%figure; raw_in(2).draw()

% discard stimuli not in use
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_channel1'; 
    'stim_channel2'; 
    'out_during';
    };
raw_in=j.run(raw_in);

% remove stimless files
j = nirs.modules.RemoveStimless();
raw_in = j.run(raw_in);


% label the short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15; % check units 
raw_in = j.run(raw_in);

% check that short-channels are labeld (command might produce error if the
% field ShortSeparation does not exist
raw_in(3).probe.link.ShortSeperation == 1 % some 1s should appear

%% Signal quality check via QT-nirs [IN]

% Extract task period and evaluate quality metrics over task period only
raw_in_ext = extract_task_period(raw_in);

j = nirs.modules.QT();
j.qThreshold = 0.6;
j.sciThreshold = 0.5;  % This paper: 10.1117/1.NPh.9.1.015001 used 0.5 as threshold
j.pspThreshold = 0; %.5; % check with psp threshold = 0 
j.fCut = [0.5 3]; % corresponds to HR between 30 and 180
ScansQuality_in = j.run(raw_in_ext);

%ScansQuality.probe.defaultdrawfcn = '10-20';

% View group results
ScansQuality_in.drawGroup('sq');
ScansQuality_in.drawGroup('sci');
ScansQuality_in.drawGroup('psp');
ScansQuality_in.drawGroup('bar');
ScansQuality_in.drawGroup('sqmask');

% View individual results 
for ii = 1:numel(ScansQuality_in)
   %figure(ii),
   ScansQuality_in(ii).drawGroup('sq');
   ScansQuality_in(ii).drawGroup('sci'); 
   ScansQuality_in(ii).drawGroup('psp');
end

%ScansQuality.qMats
%ScansQuality.probe

% Identify how many channels are bad (percentage) and where they are
% located

% Declare variables that will hold overall survival rate
ChTotALL = 0;
badChTotALL = 0;
ChLongALL = 0;
badChLongALL = 0;
ChssALL = 0;
badChssALL = 0;


fprintf('Survival rate:\n');
for i=1:1:numel(ScansQuality_in)
    % find total number of channels
    Ch = raw_in_ext(i).probe.link;
    ChTot = height(Ch)/2; % divide by two to account for hbo/hbr
    
    % find bad channels
    idxBadCh = find(ScansQuality_in(i).qMats.MeasListAct==0);
    badChTot = length(idxBadCh)/2;
    
    % find long channels
    idxChLong = find(raw_in_ext(i).probe.link.ShortSeperation == 0);
    ChLong = length(idxChLong)/2;
    
    % find bad long channels
    idxBadChLong = find((ScansQuality_in(i).qMats.MeasListAct==0) & (raw_in_ext(i).probe.link.ShortSeperation == 0));
    badChLong = length(idxBadChLong)/2;
    
    % find short channels
    idxChss = find(raw_in_ext(i).probe.link.ShortSeperation == 1);
    Chss = length(idxChss)/2;
    
    % find bad short channels
    idxBadChss = find((ScansQuality_in(i).qMats.MeasListAct==0) & (raw_in_ext(i).probe.link.ShortSeperation == 1));
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

%% Preprocessing IN during

% Fix NaNs (necessary)
j = nirs.modules.FixNaNs();
raw_out = j.run(raw_in);

% Resample to 5 HZ
j = nirs.modules.Resample();
j.Fs = 5;
raw_out = j.run(raw_out);

% Convert to optical density
j = nirs.modules.OpticalDensity();
od = j.run(raw_out);

% Wavelet fileter
j=nirs.modules.WaveletFilter(); 
odW=j.run(od);

% TDDR filter 
j = nirs.modules.TDDR(); % should be applied to od values, its preferable to not downsample before applying this algorthm 
TDDR = j.run(od);
TDDRW = j.run(odW);

% TDDR with pca option 
j = nirs.modules.TDDR();
j.usePCA = 1;
TDDRpca = j.run(od);

% Convert to hemoglobin concentration changes
j = nirs.modules.BeerLambertLaw();
hb = j.run(od);
hbW =j.run(odW);
hbTDDR = j.run(TDDR);
hbTDDRW = j.run(TDDRW);
hbTDDRpca = j.run(TDDRpca);


j = eeg.modules.BandPassFilter(); % Use high-pass filter to remove drift in signal
    j.lowpass = [];
    j.highpass = .01;
    j.do_downsample = false;
hbH = j.run(od);
hbWH =j.run(odW);
hbTDDRH = j.run(TDDR);
hbTDDRWH = j.run(TDDRW);
hbTDDRpcaH = j.run(TDDRpca);

% Label short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15;

hb = j.run(hb);
hbW = j.run(hbW);
hbTDDR = j.run(hbTDDR);
hbTDDRW = j.run(hbTDDRW);
hbTDDRpca = j.run(hbTDDRpca);


% Channel pruning based on QTNIRS results (set all bad channels to nan values)

hb_pruned = hbTDDRW;
% Channel pruning based on QT-NIRS results:
% Set bad channels to zero (both long and short) so the ShortDistanceFilter
% will work later on. 
sChannel = find(hb_pruned(1).probe.link.ShortSeperation == 1);
for i=1:numel(hb_pruned)
   idxBadCh = find(ScansQuality_in(i).qMats.MeasListAct==0);
   fprintf('Scan:%i #BadChannels:%i\n',i,length(idxBadCh)/2);
   hb_pruned(i).data(:,idxBadCh) = 0; % or nan % test both
   hb_pruned(i).data(:,sChannel(isnan(hb_pruned(i).data(1,sChannel)))) = 0; 
end 


% Use short-channels as pre-filter for fc
j = advanced.nirs.modules.ShortDistanceFilter();
hbTDDRSSfilt_in = j.run(hbTDDR);


%% Add demographics field based on condition
raw = [raw_qr raw_in raw_out_comb]
for i =1:numel(raw)
    % rename experiment variable 
    newExpStr = erase(raw(i).demographics.experiment, '_during');
    newExpStr = strip(newExpStr,'left','_'); % remove leading underscores from string
    raw(i).demographics.experiment = newExpStr;
    
  
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
end


%% 3. loading OUT during
% and Combining OUT during data files

rootdir = 'Z:\NIRS_Projects\SpecializedFoundation\BREXFIT_S3\Data\NIRS';
out_dirlist = dir(fullfile(rootdir, '**\out_during\2022*')); % get folders of outdoor during data

% First, find which files belong together
% Then, combine those two files
% Put the combined file in a new data set 

raw_out_comb = [];
combined_files = {};
kk = 1;
for i=1:1:length(out_dirlist) % iterate through dirlist 
    exp = '\d{4}'; %matches four consecutive digits.

    participant_i = char(regexp(out_dirlist(i).folder, exp, 'match'));
    folder_i = out_dirlist(i).name;
    path_i = append(out_dirlist(i).folder, '\', out_dirlist(i).name); % find full path of folder

    for j=2:1:length(out_dirlist) % iterate through dirlist again

        participant_j = char(regexp(out_dirlist(j).folder, exp, 'match'));
        folder_j = out_dirlist(j).name;
        path_j = append(out_dirlist(j).folder, '\', out_dirlist(j).name);

        if (strcmp(participant_j,participant_i) ==1) & (strcmp(folder_j,folder_i) == 0)
            fprintf('equal participant nr, different folders, --combine these!!\n');
            % find out which one is r and not
            TF = contains(folder_j,'r'); % returns 1 if folder_j contains r
            if TF == 1
                recovered_path = path_j;
                online_path = path_i;
            else
                online_path = path_j;
                recovered_path = path_i;
            end
            % combine if they have not already been combined
            if ismember(participant_i, combined_files)
                %disp(ismember(participant_i, combined_files));
                %disp('these files have been combined')
                continue
            else  
                raw_combined = combine_outdoor_files(online_path, recovered_path);
                raw_out_comb = [raw_out_comb raw_combined];
                combined_files{1, kk} = participant_i;
                kk = kk +1;
            end
        else
            continue
        end
    end
end
% This works 
%raw_out_comb.draw()
% check that probe looks alright as well
figure;raw_out_comb(1).probe.draw()

%% correct stimuli name and duration [out during]

% last timestamp minus 5 minutes and then fifteen minutes before that. So
% 20 minutes from the last timestamp, and 15 min duration. Except if the
% duration is shorter than expected (25 min), due to it ending prematurely.
% create a new stimuli marker if there isn't one

% it might be easier to erase all markers and just add new ones
%raw_out_comb_copy = raw_out_comb

for i= 1:1:length(raw_out_comb)
    last_timestamp = raw_out_comb(i).time(end);
    if last_timestamp > 20*60
        % subtract 20 minutes from last timestamp and set that as onset
        onset = last_timestamp - 20*60;
        dur = 60*15;
        name = 'out_during';
        amp = 1;
        out_during = nirs.design.StimulusEvents(name, onset, dur, amp);
        % add stimuli event to nirs file
        raw_out_comb(i).stimulus('out_during') = out_during;
        disp('recording longer than 20 min ')
    else
        % if the recording is less than 20 minutes long, 
        % skipp it 
        disp('skipping short session')
        continue
        % TODO: Decide how long these stimuli markers should be
    end
end
% check that it worked
figure; raw_out_comb(3).draw()
raw_out_comb.draw()

% discard stimuli not in use
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_channel1';
    'stim_channel2';
    'stim_aux1'; 
    'stim_aux2',
    };
raw_out_comb=j.run(raw_out_comb);

% remove stimless files
j = nirs.modules.RemoveStimless();
raw_out_comb = j.run(raw_out_comb);


% label the short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15; % check units 
raw_out_comb = j.run(raw_out_comb);

% check that short-channels are labeld (command might produce error if the
% field ShortSeparation does not exist
raw_out_comb(3).probe.link.ShortSeperation == 1 % some 1s should appear
%% Notes from QC NIRS excel sheet
% 3023: Start file is there but no initial trigger was placed.

%% Signal quality check via QT-nirs [OUT]
% Perform individual signal quality check for qr, in and out. 

% Extract task period and evaluate quality metrics over task period only
raw_out_ext = extract_task_period(raw_out_comb);

% Quality Check
j = nirs.modules.QT();
j.qThreshold = 0.6; 
j.sciThreshold = 0.5; % scalp coupling index. This paper: 10.1117/1.NPh.9.1.015001 used 0.5 as threshold
j.pspThreshold = 0; %.5; % check with psp threshold = 0 
j.fCut = [0.5 3]; % corresponds to HR between 30 and 180
ScansQuality_out = j.run(raw_out_ext);

%ScansQuality_out.probe.defaultdrawfcn = '10-20';

% View group results
ScansQuality_out.drawGroup('sq');
ScansQuality_out.drawGroup('sci');
ScansQuality_out.drawGroup('psp');
ScansQuality_out.drawGroup('bar');
ScansQuality_out.drawGroup('sqmask');

% View individual results 
for ii = 1:numel(ScansQuality_out)
   %figure(ii),
   ScansQuality_out(ii).drawGroup('sq');
   ScansQuality_out(ii).drawGroup('sci'); 
   ScansQuality_out(ii).drawGroup('psp');
end

%ScansQuality_out.qMats
%ScansQuality_out.probe

% Identify how many channels are bad (percentage) and where they are
% located

% Declare variables that will hold overall survival rate
ChTotALL = 0;
badChTotALL = 0;
ChLongALL = 0;
badChLongALL = 0;
ChssALL = 0;
badChssALL = 0;

fprintf('Survival rate:\n');
for i=1:1:numel(ScansQuality_out)
    % find total number of channels
    Ch = raw_out_ext(i).probe.link;
    ChTot = height(Ch)/2; % divide by two to account for hbo/hbr
    
    % find bad channels
    idxBadCh = find(ScansQuality_out(i).qMats.MeasListAct==0);
    badChTot = length(idxBadCh)/2;
    
    % find long channels
    idxChLong = find(raw_out_ext(i).probe.link.ShortSeperation == 0);
    ChLong = length(idxChLong)/2;
    
    % find bad long channels
    idxBadChLong = find((ScansQuality_out(i).qMats.MeasListAct==0) & (raw_out_ext(i).probe.link.ShortSeperation == 0));
    badChLong = length(idxBadChLong)/2;
    
    % find short channels
    idxChss = find(raw_out_ext(i).probe.link.ShortSeperation == 1);
    Chss = length(idxChss)/2;
    
    % find bad short channels
    idxBadChss = find((ScansQuality_out(i).qMats.MeasListAct==0) & (raw_out_ext(i).probe.link.ShortSeperation == 1));
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


%% Preproccesing [out during]
% spline interpolation and wavelet transform for motion arteafct corection before fc analysis as
% reccomended by gio, aaron and eric

% compare several pipelines
% - spline interpolation + wavelet
% tddr 
% using short channels as a prefilter for fc

% Basic pipelines
%  {1×1 nirs.modules.ImportData    }
%     {1×1 nirs.modules.RemoveStimless}
%     {1×1 nirs.modules.FixNaNs       }
%     {1×1 nirs.modules.Resample      }
%     {1×1 nirs.modules.OpticalDensity}
%     {1×1 nirs.modules.BeerLambertLaw}
%     {1×1 nirs.modules.TrimBaseline  }
%     {1×1 nirs.modules.ExportData    }

% test pipeline on copy:
%raw_out = raw_out_comb;
raw_out = raw
% Fix NaNs (necessary)
j = nirs.modules.FixNaNs();
raw_out = j.run(raw_out);

% Resample to 5 HZ
j = nirs.modules.Resample();
j.Fs = 5;
raw_out = j.run(raw_out);

% Convert to optical density
j = nirs.modules.OpticalDensity();
od = j.run(raw_out);

% Wavelet fileter
j=nirs.modules.WaveletFilter(); 
odW=j.run(od);

% TDDR filter 
j = nirs.modules.TDDR(); % should be applied to od values, its preferable to not downsample before applying this algorthm 
TDDR = j.run(od);
TDDRW = j.run(odW);

% TDDR with pca option 
j = nirs.modules.TDDR();
j.usePCA = 1;
TDDRpca = j.run(od);

% Convert to hemoglobin concentration changes
j = nirs.modules.BeerLambertLaw();
hb = j.run(od);
hbW =j.run(odW);
hbTDDR = j.run(TDDR);
hbTDDRW = j.run(TDDRW);
hbTDDRpca = j.run(TDDRpca);


j = eeg.modules.BandPassFilter(); % Use high-pass filter to remove drift in signal
    j.lowpass = [];
    j.highpass = .01;
    j.do_downsample = false;
hbH = j.run(od);
hbWH =j.run(odW);
hbTDDRH = j.run(TDDR);
hbTDDRWH = j.run(TDDRW);
hbTDDRpcaH = j.run(TDDRpca);

% Label short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15;

hb = j.run(hb);
hbW = j.run(hbW);
hbTDDR = j.run(hbTDDR);
hbTDDRW = j.run(hbTDDRW);
hbTDDRpca = j.run(hbTDDRpca);




% Channel pruning based on QTNIRS results (set all bad channels to nan values)

hb_pruned = hbTDDRW;
% Channel pruning based on QT-NIRS results:
% Set bad channels to zero (both long and short) so the ShortDistanceFilter
% will work later on. 
sChannel = find(hb_pruned(1).probe.link.ShortSeperation == 1);
for i=1:numel(hb_pruned)
   idxBadCh = find(ScansQuality(i).qMats.MeasListAct==0);
   fprintf('Scan:%i #BadChannels:%i\n',i,length(idxBadCh)/2);
   hb_pruned(i).data(:,idxBadCh) = 0; % or nan % test both
   hb_pruned(i).data(:,sChannel(isnan(hb_pruned(i).data(1,sChannel)))) = 0; 
end 


% Use short-channels as pre-filter for fc
j = advanced.nirs.modules.ShortDistanceFilter();
hbTDDRSSfilt = j.run(hbTDDR);



%% Visual comparison of different preprocessing pipelines

for i=1:numel(raw_out)
    figure;
    t = tiledlayout(3,4);

    nexttile; hb(i).draw(); title('hb (No filter)');    
    nexttile; hbH(i).draw(); title('hb highpass filter');    
    nexttile; hbTDDR(i).draw(); title('hb TDDR');    
    nexttile; hbTDDRH(i).draw(); title('hb TDDR + highpass');    
    nexttile; hbTDDRW(i).draw(); title('hb TDDR + Wavelet');
    nexttile; hbTDDRWH(i).draw(); title('hb TDDR + Wavelet + highpass');
    nexttile; hbTDDRpca(i).draw(); title('hb TDDR pca option')
    nexttile; hbTDDRpcaH(i).draw(); title('hb TDDR pca option + highpass')
    nexttile; hbWH(i).draw(); title('hb Wavelet + highpass');    
    nexttile; hbW(i).draw(); title('hb Wavelet');

    title(t, 'Hb for subject:', char(hb(i).demographics.subject))
end

%% Concatenate qr, in and out to one data set 

raw_during = [raw_qr_ext raw_in_ext raw_out_ext ];
duringDemoTable = nirs.createDemographicsTable(raw_during);
disp(duringDemoTable)   

%% Quality Check (QT-NIRS) of entire data set 

j = nirs.modules.QT();
j.qThreshold = 0.6; 
j.sciThreshold = 0.5; % scalp coupling index. This paper: 10.1117/1.NPh.9.1.015001 used 0.5 as threshold
j.pspThreshold = 0.05; %.5; % check with psp threshold = 0 
j.fCut = [0.5 3]; % corresponds to HR between 30 and 180
ScansQuality_during = j.run(raw_during);

% View group results
ScansQuality_during.drawGroup('sq');
ScansQuality_during.drawGroup('sci');
ScansQuality_during.drawGroup('psp');
ScansQuality_during.drawGroup('bar');
ScansQuality_during.drawGroup('sqmask');

% Declare variables that will hold overall survival rate
ChTotALL = 0;
badChTotALL = 0;
ChLongALL = 0;
badChLongALL = 0;
ChssALL = 0;
badChssALL = 0;

fprintf('Survival rate:\n');
for i=1:1:numel(ScansQuality_during)
    % find total number of channels
    Ch = raw_during(i).probe.link;
    ChTot = height(Ch)/2; % divide by two to account for hbo/hbr
    
    % find bad channels
    idxBadCh = find(ScansQuality_during(i).qMats.MeasListAct==0);
    badChTot = length(idxBadCh)/2;   
    % find long channels
    idxChLong = find(raw_during(i).probe.link.ShortSeperation == 0);
    ChLong = length(idxChLong)/2;
    % find bad long channels
    idxBadChLong = find((ScansQuality_during(i).qMats.MeasListAct==0) & (raw_during(i).probe.link.ShortSeperation == 0));
    badChLong = length(idxBadChLong)/2;    
    % find short channels
    idxChss = find(raw_during(i).probe.link.ShortSeperation == 1);
    Chss = length(idxChss)/2;   
    % find bad short channels
    idxBadChss = find((ScansQuality_during(i).qMats.MeasListAct==0) & (raw_during(i).probe.link.ShortSeperation == 1));
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
% Survival rate:
% Scan:1 #ChannelsTotal: 69.32 percent #ChannelsLong: 65.28 percent #ChannelsShort: 87.50 percent
% Scan:2 #ChannelsTotal: 61.36 percent #ChannelsLong: 58.33 percent #ChannelsShort: 75.00 percent
% Scan:3 #ChannelsTotal: 54.55 percent #ChannelsLong: 61.11 percent #ChannelsShort: 25.00 percent
% Scan:4 #ChannelsTotal: 31.82 percent #ChannelsLong: 36.11 percent #ChannelsShort: 12.50 percent
% Scan:5 #ChannelsTotal: 63.64 percent #ChannelsLong: 66.67 percent #ChannelsShort: 50.00 percent
% Scan:6 #ChannelsTotal: 56.82 percent #ChannelsLong: 62.50 percent #ChannelsShort: 31.25 percent
% Scan:7 #ChannelsTotal: 52.27 percent #ChannelsLong: 52.78 percent #ChannelsShort: 50.00 percent
% Scan:8 #ChannelsTotal: 7.95 percent #ChannelsLong: 8.33 percent #ChannelsShort: 6.25 percent
% Scan:9 #ChannelsTotal: 23.86 percent #ChannelsLong: 26.39 percent #ChannelsShort: 12.50 percent
% Scan:10 #ChannelsTotal: 75.00 percent #ChannelsLong: 79.17 percent #ChannelsShort: 56.25 percent
% Scan:11 #ChannelsTotal: 81.82 percent #ChannelsLong: 84.72 percent #ChannelsShort: 68.75 percent
% Scan:12 #ChannelsTotal: 89.77 percent #ChannelsLong: 94.44 percent #ChannelsShort: 68.75 percent
% Scan:13 #ChannelsTotal: 90.91 percent #ChannelsLong: 91.67 percent #ChannelsShort: 87.50 percent
% Scan:14 #ChannelsTotal: 17.05 percent #ChannelsLong: 18.06 percent #ChannelsShort: 12.50 percent
% Scan:15 #ChannelsTotal: 55.68 percent #ChannelsLong: 56.94 percent #ChannelsShort: 50.00 percent
% Scan:16 #ChannelsTotal: 71.59 percent #ChannelsLong: 68.06 percent #ChannelsShort: 87.50 percent
% Scan:17 #ChannelsTotal: 31.82 percent #ChannelsLong: 31.94 percent #ChannelsShort: 31.25 percent
% Scan:18 #ChannelsTotal: 30.68 percent #ChannelsLong: 33.33 percent #ChannelsShort: 18.75 percent
% Scan:19 #ChannelsTotal: 43.18 percent #ChannelsLong: 41.67 percent #ChannelsShort: 50.00 percent
% Scan:20 #ChannelsTotal: 20.45 percent #ChannelsLong: 11.11 percent #ChannelsShort: 62.50 percent
% Scan:21 #ChannelsTotal: 45.45 percent #ChannelsLong: 44.44 percent #ChannelsShort: 50.00 percent
% Scan:22 #ChannelsTotal: 92.05 percent #ChannelsLong: 95.83 percent #ChannelsShort: 75.00 percent
% Scan:23 #ChannelsTotal: 35.23 percent #ChannelsLong: 38.89 percent #ChannelsShort: 18.75 percent
% Scan:24 #ChannelsTotal: 75.00 percent #ChannelsLong: 76.39 percent #ChannelsShort: 68.75 percent
% Scan:25 #ChannelsTotal: 80.68 percent #ChannelsLong: 80.56 percent #ChannelsShort: 81.25 percent
% Scan:26 #ChannelsTotal: 90.91 percent #ChannelsLong: 93.06 percent #ChannelsShort: 81.25 percent
% Scan:27 #ChannelsTotal: 36.54 percent #ChannelsLong: 36.36 percent #ChannelsShort: 37.50 percent
% Scan:28 #ChannelsTotal: 30.77 percent #ChannelsLong: 29.55 percent #ChannelsShort: 37.50 percent
% Scan:29 #ChannelsTotal: 11.54 percent #ChannelsLong: 9.09 percent #ChannelsShort: 25.00 percent
% Scan:30 #ChannelsTotal: 25.00 percent #ChannelsLong: 27.27 percent #ChannelsShort: 12.50 percent
% Scan:31 #ChannelsTotal: 30.77 percent #ChannelsLong: 34.09 percent #ChannelsShort: 12.50 percent
% Scan:32 #ChannelsTotal: 9.62 percent #ChannelsLong: 6.82 percent #ChannelsShort: 25.00 percent
% Scan:33 #ChannelsTotal: 32.69 percent #ChannelsLong: 34.09 percent #ChannelsShort: 25.00 percent
% Scan:34 #ChannelsTotal: 13.46 percent #ChannelsLong: 13.64 percent #ChannelsShort: 12.50 percent
% Scan:35 #ChannelsTotal: 61.54 percent #ChannelsLong: 59.09 percent #ChannelsShort: 75.00 percent
% Scan:36 #ChannelsTotal: 32.69 percent #ChannelsLong: 29.55 percent #ChannelsShort: 50.00 percent
% Overall survival rate: #ChannelsTotal: 50.68 percent #ChannelsLong: 51.34 percent #ChannelsShort: 47.58 percent
%% Preprocessing [entire data set]
% work on a copy
raw_during_ppp = raw_during;

% Fix NaNs (necessary)
j = nirs.modules.FixNaNs();
raw_during_ppp = j.run(raw_during_ppp);

% Resample to 5 HZ
j = nirs.modules.Resample();
j.Fs = 5;
raw_during_ppp = j.run(raw_during_ppp);

% Convert to optical density
j = nirs.modules.OpticalDensity();
od = j.run(raw_during_ppp);

% Wavelet fileter
j=nirs.modules.WaveletFilter(); 
odW=j.run(od);

% TDDR filter (with and without additional wavelet filter)
j = nirs.modules.TDDR(); % should be applied to od values, its preferable to not downsample before applying this algorthm 
TDDR = j.run(od);
WTDDR = j.run(odW);

% TDDR with pca option 
j = nirs.modules.TDDR();
j.usePCA = 1;
TDDRpca = j.run(od);

% TDDR first, then wavelet filter
j=nirs.modules.WaveletFilter(); 
TDDRW=j.run(TDDR);

% Bandpassfilter
j = eeg.modules.BandPassFilter(); % Use high-pass filter to remove drift in signal
    j.lowpass = []; % m
    j.highpass = .01;
    j.do_downsample = false;

odH = j.run(od);
odWH =j.run(odW);
TDDRH = j.run(TDDR);
TDDRWH = j.run(TDDRW);
WTDDRH = j.run(WTDDR);
TDDRpcaH = j.run(TDDRpca);


% Convert to hemoglobin concentration changes
j = nirs.modules.BeerLambertLaw();
hbWTDDR = j.run(WTDDR);

hb = j.run(od);
hbW =j.run(odW);
hbTDDR = j.run(TDDR);
hbTDDRW = j.run(TDDRW);

hbTDDRpca = j.run(TDDRpca);

hbH = j.run(odH);
hbodWH =j.run(odWH);
hbTDDRH = j.run(TDDRH);
hbTDDRWH = j.run(TDDRWH);
hbWTDDRH = j.run(WTDDRH);
hbTDDRpcaH = j.run(TDDRpcaH);



% Label short-channels
j=nirs.modules.LabelShortSeperation();
j.max_distance=15;

hb = j.run(hb);
hbW = j.run(hbW);
hbTDDR = j.run(hbTDDR);
hbTDDRW = j.run(hbTDDRW);
%hbWTDDR = j.run(hbWTDDR);
hbTDDRpca = j.run(hbTDDRpca);

hbH = j.run(hbH);
hbodWH =j.run(hbodWH);
hbTDDRH = j.run(hbTDDRH);
hbTDDRWH = j.run(hbTDDRWH);
hbWTDDRH = j.run(hbWTDDRH);
hbTDDRpcaH = j.run(hbTDDRpcaH);


% Channel pruning based on QT-NIRS results:
hb_pruned = prune_channels(hb,ScansQuality_during);
hbW_pruned = prune_channels(hbW,ScansQuality_during);
hbTDDR_pruned = prune_channels(hbTDDR,ScansQuality_during);
hbTDDRW_pruned = prune_channels(hbTDDRW,ScansQuality_during);

hbWTDDR_pruned = prune_channels_nan(hbWTDDR,ScansQuality_during);

hbTDDRpca_pruned = prune_channels(hbTDDRpca,ScansQuality_during);
hbH_pruned = prune_channels(hbH,ScansQuality_during);
hbodWH_pruned =prune_channels(hbodWH,ScansQuality_during);
hbTDDRH_pruned =prune_channels(hbTDDRH,ScansQuality_during);
hbTDDRWH_pruned = prune_channels(hbTDDRWH,ScansQuality_during);
hbWTDDRH_pruned = prune_channels(hbWTDDRH,ScansQuality_during);
hbTDDRpcaH_pruned = prune_channels(hbTDDRpcaH,ScansQuality_during);


% Use short-channels as pre-filter for fc
j = advanced.nirs.modules.ShortDistanceFilter();
hb_prunedSSfilt = j.run(hb_pruned);
hbWTDDRH_prunedSSfilt = j.run(hbWTDDRH_pruned);


% Visual comparison of different preprocessing pipelines

for i=1:numel(raw_during)
    figure;
    t = tiledlayout(3,4);

    nexttile; hb_pruned(i).draw(); title('(No filter)');
    nexttile; hb_prunedSSfilt(i).draw(); title('ShortChannelFilter');
    nexttile; hbH_pruned(i).draw(); title('highpass');   

    nexttile; hbTDDR_pruned(i).draw(); title('TDDR');    
    nexttile; hbTDDRH_pruned(i).draw(); title('TDDR + highpass');    
    nexttile; hbTDDRW_pruned(i).draw(); title('TDDR + Wavelet');
    nexttile; hbWTDDR_pruned(i).draw(); title('Wavelet + TDDR');
    nexttile; hbTDDRpca_pruned(i).draw(); title('TDDR pca option')
    nexttile; hbTDDRpcaH_pruned(i).draw(); title('TDDR pca option + highpass')
    nexttile; hbodWH_pruned(i).draw(); title('Wavelet + highpass');    
    nexttile; hbW_pruned(i).draw(); title('Wavelet');
    nexttile; hbWTDDRH_prunedSSfilt(i).draw(); title('Wavelet + TDDR + ShortChannelFilter');

    title(t, 'Hb for subject:', char(hb_pruned(i).demographics.subject), char(hb_pruned(i).demographics.load))
end


for i=1:6:numel(raw_during)
    figure;
    t = tiledlayout(3,1);

    nexttile; hb_pruned(i).draw(); title('(No filter)');
    nexttile; hb_prunedSSfilt(i).draw(); title('ShortChannelFilter');
    nexttile; hbWTDDRH_prunedSSfilt(i).draw(); title('Wavelet + TDDR + ShortChannelFilter');
    string  = append('Hb for subject: ', char(hb_pruned(i).demographics.subject), ' ', char(hb_pruned(i).demographics.load));
    title(t, string);
end


for i=1:6:numel(raw)
    figure;
    t = tiledlayout(3,2);

    nexttile; hb(i).draw(); title('(No filter)');
    nexttile; hbW(i).draw(); title('wavelet');
    nexttile; hbTDDR(i).draw(); title('TDDR');
    nexttile; hbTDDRW(i).draw(); title('TDDR + Wavelet');
    nexttile; hbTDDRpca(i).draw(); title('TDDR pca option')
    string  = append('Hb for subject: ', char(hb(i).demographics.subject), ' ', char(hb(i).demographics.load));
    title(t, string);
end


for i=1:numel(raw)
    if hb.dem
    figure;
    t = tiledlayout(3,2);

    nexttile; hb(i).draw(); title('(No filter)');
    nexttile; hbW(i).draw(); title('wavelet');
    nexttile; hbTDDR(i).draw(); title('TDDR');
    nexttile; hbTDDRW(i).draw(); title('TDDR + Wavelet');
    nexttile; hbTDDRpca(i).draw(); title('TDDR pca option')
    string  = append('Hb for subject: ', char(hb(i).demographics.subject), ' ', char(hb(i).demographics.load));
    title(t, string);
end
%%


% discard stimuli not in use
j=nirs.modules.DiscardStims;
j.listOfStims={
    'stim_channel1';
    'stim_channel2';
    'stim_aux1'; 
    'stim_aux2',
    };
hbWTDDRH_prunedSSfilt=j.run(hbWTDDRH_prunedSSfilt);
%% 1st (Subject) level connectivity model [entire data set ]


j = nirs.modules.Connectivity(); 
j.corrfcn %leave default
%j.corrfcn  = @(data)nirs.sFC.ar_corr(data,'32xFs',true)
j.divide_events = 1; % fc during cycling
ConnStats = j.run(hbWTDDR_pruned);

% ConnStatsQR = j.run(conQR);
% ConnStatsIN = j.run(conIN);
% ConnStatsOUT = j.run(conOUT);

head(ConnStats(1).table(),20)

%% 2nd (Group) Level connectivity
j = nirs.modules.MixedEffectsConnectivity()
GroupConnStats = j.run(ConnStats)

conQR = [];
conIN = [];
conOUT  = [];
for i = 1:numel(ConnStats)
    if strcmp(ConnStats(i).demographics.load, 'qr')
        conQR = [conQR ConnStats(i)];
    elseif strcmp(ConnStats(i).demographics.load, 'in')
        conIN = [conIN ConnStats(i)];
    elseif strcmp(ConnStats(i).demographics.load, 'out')
        conOUT = [conOUT ConnStats(i)];
    end
end

% breaks when ran on pruned data (set to 0, try with nan)

j = nirs.modules.MixedEffectsConnectivity()
GroupConnStatsQR = j.run(conQR) 

head(GroupConnStatsQR.table())
GroupConnStatsQR.draw('t',[-5 5],'qvalue<0.05');

j = nirs.modules.MixedEffectsConnectivity()
GroupConnStatsIN = j.run(conIN) 

j = nirs.modules.MixedEffectsConnectivity()
GroupConnStatsOUT = j.run(conOUT) 

 % ##### SECOND ATTEMPT 
j = nirs.modules.Connectivity()
j.corrfcn =@(data)nirs.sFC.ar_wcoher(data,'4x',[.05 .2],'morl',true); % Whitened Wavelet coherence
j.divide_events = 1; % fc during cycling
ConnStats2 = j.run(hbSS_Filt2)

j = nirs.modules.MixedEffectsConnectivity()
GroupConnStats2 = j.run(ConnStats2)

head(GroupConnStats2.table())
figure; GroupConnStats2.draw()
%('R',[-1 1],'q<0.05')
figure; GroupConnStats2.draw('t',[-5 5],'qvalue<0.2')




%% FC analysis on outdoor data set
cleanhb(3).probe.link.ShortSeperation == 1 % check that short channels are registered


% functional connectivity during cycling
j = nirs.modules.Connectivity()
j.corrfcn %leave default
    % Use one of the two below i think?
   % j.corrfcn =@(data)nirs.sFC.ar_wcoher(data,'4x',[.05 .2],'morl',true); % Whitened Wavelet coherence
   % j.corrfcn = @(data) nirs.sFC.ar_corr(data,'32xFs',true);

j.divide_events = 1; % functional connectivity during cycling
ConnStatsDivEvArW = j.run(hbSS)
% would run on hbSS
% did not run on hbSS_filt
% need to find ot where it breaks

figure; ConnStatsDivEvArW(2).draw('R',[-1 1],'q<0.05')


j = nirs.modules.MixedEffectsConnectivity();

GroupConnStats = j.run(ConnStatsDivEvArW);
%%
% Questions:
% - Should we define ROIs for fc analysis (instead of channels)

%% single participant

% HRF?
HRF = GroupStats.HRF
nirs.viz.plot2D(HRF)
figure; HRF.draw()

% functional connectivity
% Question: How to include short channels in the fc analysis?

%nirs.modules.connectivity_GLM

j = nirs.modules.Connectivity()
j.corrfcn % use default  @(data)nirs.sFC.ar_corr(data,'4xFs',true)

ConnStats = j.run(hbResampled)


% You can draw the connectiivty maps by:
figure; ConnStats.draw('R',[-1 1],'p<0.05')

figure; ConnStats.draw('R',[-1 1],'q<0.05')

figure; ConnStats.draw('R',[],'q<0.05')

% functional connectivity during cycling
j = nirs.modules.Connectivity()
j.corrfcn % use default  @(data)nirs.sFC.ar_corr(data,'4xFs',true)
j.divide_events = 1 % functional connectivity during cycling
ConnStatsDivEv = j.run(hbResampled)
figure; ConnStatsDivEv.draw('R',[-1 1],'q<0.05')

% functional connectivity during cycling
j = nirs.modules.Connectivity()
j.corrfcn =@(data)nirs.sFC.ar_wcoher(data,'4x',[.05 .2],'morl',true); % Whitened Wavelet coherence
j.divide_events = 1 % functional connectivity during cycling
ConnStatsDivEvArW = j.run(hbResampled)
figure; ConnStatsDivEvArW.draw('R',[-1 1],'q<0.05')

figure; ConnStatsDivEvArW.draw('R',[-1 1],'q<0.005')
figure; ConnStatsDivEvArW.draw('t',[-5 5],'q<0.005')


figure; ConnStatsDivEvArW.draw('Z',[],'q<0.05')

% Cassies fc script:
% https://osf.io/tj6xr?view_only=0c977c2760054e11a888576feb55745a
% nirs-toolbox fc demo:
% https://github.com/huppertt/nirs-toolbox/blob/master/demos/fnirs_connectivity_demo.m

% fNIRS Analysis Club Facebook
% Hi, Any suggestions for task-based functional connectivity analysis?
% Which toolbox or stat approach do you recommend? % 
% Hendrik Santosa 
% You can use connectivity module (see
% fNIRS_connectivity_demo.m for demo). If you have multiple conditions, you
% will get R, Z etc for every conditions in your stats results (by default
% will be set to "Rest"). Please change the setting in the connectivity
% module: divided_events (parse to multiple conditions) and
% min_event_duration.
% Ali Rahimpour
% Thank you, Hendrik! Yes, I have
% been trying to make changes in NIRS toolbox based on my need of task
% based approach. So you meant that it still would be fine/valid to use
% stat analysis (for example AR-W) for task- based approach as the way you
% did in resting state data demo. Correct? 
% Hendrik Santosa 
% Yes correct, you can do like that

% We can also convert the connectivity object into a weighted graph object
Graph=ConnStats.graph('Z:hbo','p<0.005');

Graph.pagerank

PgRankGraph = Graph.pagerank
PgRankGraph.nodeInfo

figure; Graph.pagerank.draw;

% this breaks: "Number of observations N must be >= 2." so probably because
% we have only one participant
% j = nirs.modules.MixedEffectsConnectivity();
% GroupConnStats = j.run(ConnStats);

%% correlation matrix
stats = GroupConnStats2.R(:,:,1);
imagesc(stats) 
colorbar
%title = 
plotTitle = 'Prefrontal Cortex Resting State Connectivity';
cLabel = 'R-values';
cAxis = [0 1];
color = 'jet';
plotTitle, cLabel, cAxis, color);
plot(stats, plotTitle, cLabel, cAxis, color)
correlationMatrixLabeling(stats, plotTitle, cLabel, cAxis, color);

%% correlation matrix
stats = GroupConnStats2.Z(:,:,1);
zimg = imagesc(stats)
plotTitle = 'Prefrontal Cortex Resting State Connectivity';
cLabel = 'Z-values';
cAxis = [0 1];
color = 'jet';
correlationMatrixLabeling(stats, plotTitle, cLabel, cAxis, color);
