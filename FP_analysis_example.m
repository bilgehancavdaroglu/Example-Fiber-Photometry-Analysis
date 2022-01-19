% Fiber Photometry Analysis example
%
% First extract the data from the appropriate folder, creating a structure
% FP.SubjectID  - > Subject IDs, Timestamp can have BatchID
% FP.Experiment  - > Name of the Experiment.
% FP.RawData  - > exported from csv files as a table

tic;
clear all
cd('C:\Users\Bilge\Google Drive\Ito Lab\Rat Experiments\FP\Nisma\Appetitive_Arm MPFC');
t = natsortfiles(dir); %natsort makes sure 10 comes after 9 or 100 comes after 99 etc.
n = 0;
for iFile = 1:length(t) %For each file in the folder
    if ~t(iFile).isdir %Ignore the directory files ('.' and '..')
        n = n + 1;
        
        C = strsplit(t(iFile).name);%Split the filename to the components
        
        FP.SubjectID(n) = C(1);%This info depends on the naming of the file
        FP.Experiment(n) = C(2);%Name of the experiment
        FP.Region{n} = C{1, 3}(1:end - 4);%Brain region
        
        if  FP.Experiment(n)~="Timestamps" %Only the csv files
            FP.RawData{n} = readtable(t(iFile).name);%reads csv
        else
            FP.Timestamps = readtable(t(iFile).name);%reads the timestamps
        end
        fclose('all');%make sure files are not open in matlab
    end
end

clearvars  - except FP
toc;

%% Applying normalization methods
tic
for iSubj = 1:length(FP.SubjectID) - 1% - 1 for the timestamp
    
    raw_reference = FP.RawData{iSubj}{1:end - 1, 3}';%Make sure it is the correct column, depends on the doric rig setup alternatively you can index with the variable name i.e. AnalogIn__Ch_1_1
    raw_signal = FP.RawData{iSubj}{1:end - 1, 2}';%Make sure it is the correct column, depends on the doric rig setup
    
    % Martinova Method
    % Default option of this method ignores the first 200 lines, can be
    % adjusted, also smoothing can be changed
    FP.NormData.Martinova{iSubj}(:, 1) = FP.RawData{iSubj}{201:end, 1};%This method ignores the first 200 rows
    FP.NormData.Martinova{iSubj}(:, 2) = get_zdFF(raw_reference, raw_signal)';%Actual function (get_zdFF) that applies the martinova normalization
    
    %'DeltaFF Method
    % Data before the window won't be reliable. removed the first 200 data
    % points to compare it to the Martinova method
    binSize = mean(diff(FP.RawData{iSubj}{:, 1}));%Average binSize of the
    %data. Should be around 0.82ms, but depends on the rig setup
    window = 5; %window for the deltaF / F calculations, in seconds
    
    FP.NormData.Doric{iSubj}(:, 1) = FP.RawData{iSubj}{201:end, 1};%This method ignores the first 200 rows to make it consistent with Martinova
    FP.NormData.Doric{iSubj}(:, 2) = zscore(deltaFF3(raw_signal(200:end)', binSize, window)...
        - deltaFF3(raw_reference(200:end)', binSize, window)); %z scored
end
toc

%% Calculating the average Amplitude and AUC before and after Appetitive arm entry
%
% Use timestamp information to select events of interest during the experiment
%
% Use the normalized data for the specific analysis you want. Most likely
% candidates are using the timestamps to look at specific events, 
% calculating average amplitude and / or AUC.
% Filter the table for the relevant information.

rows = (FP.Timestamps.Arm == "Appetitive" & FP.Timestamps.Condition == "In");
tempTable = FP.Timestamps(rows, :);%Only Entries for the Appetitive Arm
range = 2; % We are interested in a 4s window with the event in the middle

for iSubj = 1:length(FP.SubjectID) - 1%    
    % timestamps for each entry / trial for the specific subject
    t = tempTable.Timestamp((tempTable.SubjID == str2double(FP.SubjectID(iSubj))), :);
    
    FP.Events(iSubj).Entry.Signal = nan(sum(~isnan(t)), 500);% Init signal
    
    for entry = 1:sum(~isnan(t)) %For each non NaN entry
        
        % Selecting the correct rows using the timestamp information, First
        % column of the NormData is the timestamps from the rig.
        ts = [find(round(FP.NormData.Martinova{iSubj}(:, 1), 2) == t(entry) - range, 1), find(round(FP.NormData.Martinova{iSubj}(:, 1), 2) == t(entry)  +  range, 1)];
        
        FP.Events(iSubj).Entry.Signal(entry, 1:(ts(2) - ts(1)  +  1)) = FP.NormData.Martinova{iSubj}(ts(1):ts(2), 2)';
        
        FP.Events(iSubj).Entry.AMP(entry, 1) = mean( FP.Events(iSubj).Entry.Signal(entry, 1:round((ts(2) - ts(1)) / 2)), 'omitnan');%Average Amplitude range seconds before
        FP.Events(iSubj).Entry.AMP(entry, 2) = mean( FP.Events(iSubj).Entry.Signal(entry, (round((ts(2) - ts(1)) / 2)) + 1:round(ts(2) - ts(1))), 'omitnan');%Average Amplitude range seconds after
        
        FP.Events(iSubj).Entry.AUC(entry, 1) = trapz( FP.Events(iSubj).Entry.Signal(entry, 1:round((ts(2) - ts(1)) / 2)));%Area under the curve range seconds before
        FP.Events(iSubj).Entry.AUC(entry, 2) = trapz( FP.Events(iSubj).Entry.Signal(entry, round((ts(2) - ts(1)) / 2) + 1:round(ts(2) - ts(1))));%Area under the curve range seconds after
        
    end
end






