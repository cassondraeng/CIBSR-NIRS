%% Example: filtering data based on QC NIRX excel sheet
% This script exemplifies how to filter fnirs data files, e.g., based on 
% notes from the QC-NIRS excel sheet. I added a separate script for this to
% a) remove a redundant filtering section from the nback-script, and b)
% because it might be nice to have to have the option of looking at the
% filtering separate from loading the files to better understand what's
% going on. 

% To run this script you need nirs-toolbox since the filtering assumes that
% the finirs data is stored in the toolbox's data format. 

% Last updated: 2022-12-19
% Created by: Henrikke Dybvik
% Using MATLAB R2021b

%% go though data set based on NIRS_QC.xlsx and sort/filter data

% this sorts through the data and identifies which files are the correct
% ones. It has been combined with the section above that loads, so that only
% the correct files are loaded. Therefore, it's not neccessary to run this
% section at the moment. It has been left here in case additional sorting
% becomes neccessary when more data has been collected. 

for i=1:1:numel(raw)
    participant_id = raw(i).demographics.subject;
    filedesc = raw(i).description;
    experiment = raw(i).demographics.experiment;

    if (participant_id == '3008') & startsWith(experiment, 'qr_post') & endsWith(filedesc, '005\')
        % qr post nback: two files: 004=5.2mb, 005=8.1mb -> File 005 is correct full data file
        filedesc % correct file, remove the other one 
    elseif (participant_id == '3016') & startsWith(experiment, 'qr_post') & endsWith(filedesc, '006\')
        % two files: 005=227kb, 006=8.5, File 006 is correct full data file
        filedesc % correct file, remove the other one
    elseif (participant_id == '3014') & startsWith(experiment, 'out_post') & endsWith(filedesc, '004\')
        % two files: 003=839kb, 004=8mb, File 004 is correct full data file
        filedesc
    end

end  

% go though data set based on NIRS_QC.xlsx and delete data
% this deletes files from raw, but it produces an error since raw changes
% size during iteration. Thus it's better to save the index of this file
% and remove those files after the loop is done. 

for i=1:1:numel(raw)
    participant_id = raw(i).demographics.subject;
    filedesc = raw(i).description;
    experiment = raw(i).demographics.experiment;

    if (participant_id == '3008') & startsWith(experiment, 'qr_post') & endsWith(filedesc, '004\')
        % qr post nback: two files: 004=5.2mb, 005=8.1mb -> File 005 is correct full data file
        filedesc % remove this one
        i
        raw(i) = []
    elseif (participant_id == '3016') & startsWith(experiment, 'qr_post') & endsWith(filedesc, '005\')
        % two files: 005=227kb, 006=8.5, File 006 is correct full data file
        filedesc % remove this one
        i
        raw(i) = []
    elseif (participant_id == '3014') & startsWith(experiment, 'out_post') & endsWith(filedesc, '003\')
        % two files: 003=839kb, 004=8mb, File 004 is correct full data file
        filedesc
        i
        raw(i) = []
    end
end  