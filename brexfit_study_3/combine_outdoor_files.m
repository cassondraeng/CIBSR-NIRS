function [raw_combined] = combine_outdoor_files(online_path,recovered_path)
% This function combines two fnris data files from outdoor during.
% This function takes in the path to two folders that contain outdoor
% during fnris data files. These paths should correspond to the online and 
% recovered file from the one individual participant.
    
    % combining files
    % load data into raw structure
    r1 = nirs.io.loadDirectory(online_path); % Online files is labeled 1 
    r2 = nirs.io.loadDirectory(recovered_path); % Recovered file is labeled 2
    % combine them into r3. 
    r3 = r1; % start with online file
    
    % concatenate fnirs time series data vertically
    timeseries = [r1.data; r2.data];
    % replace time series data
    r3.data= timeseries;
    % check if it worked
    %r3 % yes, theere is a missmatch between number of data points and time stamps
    
    % time is in seconds
    t1 = r1.time;
    t2 = r2.time;
    
    %size(t1)
    %size(t2)
    
    %get the last value for t1
    %t1(end,1)
    % add this value to all times in t2
    %t2 + t1(end,1)
    % this gives two rows with the same timestamp. this might pose a problem
    % later on, but has so far not
    
    t3 = [t1; t2 + t1(end,1)]; % concatenate the two timestamps
    % check dimentions
    size(t3);
    size(r3.data);
    % replace timestamps
    r3.time = t3;
    
    % visual inspection 
    %figure; r3.draw()

    raw_combined = r3;
end