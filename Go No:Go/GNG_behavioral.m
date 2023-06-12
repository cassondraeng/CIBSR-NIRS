addpath /Volumes/Projects/FragileX/FragileX_AAA/Imaging/T3/
addpath /Users/cassondraeng/Documents/MATLAB/GNG

cd /Volumes/Projects/FragileX/FragileX_AAA/Imaging/T3/

allParticipants = dir('*/*/*.mat');

avgRT = [];
avgAcc = [];
pID = cell(length(allParticipants),1);
for p = 1:length(allParticipants)

    dirPi = allParticipants(p).folder;
    filePi = allParticipants(p).name;
    behDataPi = load([dirPi, '/', filePi]);

    avgAcc = [avgAcc; behDataPi.gngHome.taskAvg];
    avgRT = [avgRT; behDataPi.gngHome.rtAvg];
    pID{p} = filePi;

end %end p loop

performanceThrIndx = find(avgAcc > 0.5); %finding thersholds

display('Keyboard End')
keyboard

%Copy avgAcc avgRT pID to Excel file