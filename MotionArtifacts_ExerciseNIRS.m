%% Script to analyze Exercise-VR -- CIBSR fNIRS Script
% What you will need before run this script: NIRS Analyzer Toolbox and customized functions
% Code from Motion Artifact and Corrections by Dr. Ted Huppert 
% Adapted by Dr. Cassie Eng 6-26-2023

clc; clear all; close all %% clear Workspace


%%

% Load the data for a condition x time design
raw = nirs.io.loadDirectory('.',{'Group','Subject','Time'});

raw %shows whats in the files 

% check variable names
nirs.createDemographicsTable(raw)

%% rename stim markers and pull blocks
raw=nirs.viz.StimUtil(raw) % change stim markers and rename the events in GUI
raw=job.run(raw);

%% visualize the data
raw.draw
raw.gui
raw(1).probe.draw
%raw(1).probe.defaultdrawfcn='?'; % all the ways you can visualize the data!
raw(1).probe.defaultdrawfcn='3D mesh (top)'; % 3D brain model
%raw(1).probe.defaultdrawfcn='10-20'; % EEG topo map
raw(1).probe=raw(1).probe.SetFiducialsVisibility(false); %removes the EEG data points
%raw(1).probe.defaultdrawfcn='3D label ball'; %I also like this one

%% Draw accelerometer 
raw(1).auxillary
acc=raw(1).auxillary('aux')
acc
acc.draw

%%
%Show short sep info:
job=nirs.modules.LabelShortSeperation(job);
job.max_distance=15;
raw = job.run(raw);

% Basic pipeline
job=nirs.modules.OpticalDensity;
job=nirs.modules.Resample(job);
job.Fs=1;
job=nirs.modules.BeerLambertLaw(job);

%% Add accelerometer data to regression model
job=nirs.modules.AddAuxRegressors(job);
job.label={'aux'};

%% Add short seperations to GLM
job=nirs.modules.GLM(job);
job.AddShortSepRegressors=true;
job = nirs.modules.KeepTypes(job); % Remove deoxy-hemoglobin
    job.types = {'hbo'}

% SUBJECT LEVEL ANALYSIS 
SubjStats=job.run(raw); % individual subject data
SubjStats 

SubjStats1(1).draw('tstat',[-5 5],'q<.05') %visual for subject 1
SubjStats1(2).draw('tstat',[-5 5],'q<.05') %visual for subject 2

%% GROUP-LEVEL ANALYSIS
% to put both files together into one structure (e.g. to compare)
job=nirs.modules.GLM(job);
job.AddShortSepRegressors=true;
job=nirs.modules.MixedEffects;
job.formula='beta ~ -1 + Group';
GroupStats=job.run(SubjStats);
GroupStats.draw('tstat', [-1 1]);
disp(GroupStats.conditions); % double check your conditions


%% OLS
job1=nirs.modules.OpticalDensity;
job1=nirs.modules.Resample(job1);
job1.Fs=1;
job1=nirs.modules.BeerLambertLaw(job1);
% job1=nirs.modules.RemoveShortSeperations(job1); BUG HERE
job1=nirs.modules.GLM(job1);
job1.type='OLS';

SubjStats1=job1.run(raw);
SubjStats1(1).draw('tstat',[-5 5],'q<.05')
SubjStats1(2).draw('tstat',[-5 5],'q<.05')

GroupStats1=job.run(SubjStats1);
GroupStats1.draw('tstat', [-5 5], 'p < 0.05');



%% AR-IRLS
job2=nirs.modules.OpticalDensity;
job2=nirs.modules.Resample(job2);
job2.Fs=1;
job2=nirs.modules.BeerLambertLaw(job2);
% job2=nirs.modules.RemoveShortSeperations(job2);BUG HERE
job2=nirs.modules.GLM(job2);

SubjStats2=job2.run(raw);
GroupStats2=job.run(SubjStats2);

%% AR-IRLS + SSreg
job3=nirs.modules.OpticalDensity;
job3=nirs.modules.Resample(job3);
job3.Fs=1;
job3=nirs.modules.BeerLambertLaw(job3);
job3=nirs.modules.GLM(job3);
job3.AddShortSepRegressors=true;

SubjStats3=job3.run(raw);
GroupStats3=job.run(SubjStats3);

%% AR-IRLS + Aux reg
job4=nirs.modules.OpticalDensity;
job4=nirs.modules.Resample(job4);
job4.Fs=1;
job4=nirs.modules.BeerLambertLaw(job4);
% job4=nirs.modules.RemoveShortSeperations(job4); BUG HERE 
job4=nirs.modules.AddAuxRegressors(job4);
job4.label={'aux'};
job4=nirs.modules.GLM(job4);

SubjStats4=job4.run(raw);
GroupStats4=job.run(SubjStats4);


%% AR-IRLS + Aux reg + SSreg
job5=nirs.modules.OpticalDensity;
job5=nirs.modules.Resample(job5);
job5.Fs=1;
job5=nirs.modules.BeerLambertLaw(job5);
job5=nirs.modules.AddAuxRegressors(job5);
job5.label={'aux'};
job5=nirs.modules.GLM(job5);
job5.AddShortSepRegressors=true;
SubjStats5=job5.run(raw);
GroupStats5=job.run(SubjStats5);

%% Comparing 5 models. 
% OC curve (receiver operating characteristic curve) shows the performance of a classification model 
% at all classification thresholds. This curve plots two parameters: True Positive Rate. False Positive Rate.
ROC1=nirs.testing.ChannelStatsROC;
ROC1.simfunc=@()nirs.testing.simData(raw(1));
ROC1.pipeline={job1 job2 job3 job4 job5};
ROC1=ROC1.run(50); 
ROC1.draw('hbo');
% Receiver Operator Curve (ROC) analysis runs 50 iterations of each model
% for Condition 1

ROC2=nirs.testing.ChannelStatsROC;
ROC2.simfunc=@()nirs.testing.simData(raw(2));
ROC2.pipeline={job1 job2 job3 job4 job5};
ROC2=ROC2.run(50);
ROC2.draw('hbo');
% Receiver Operator Curve (ROC) analysis runs 50 iterations of each model
% for Condition 2 etc. 




