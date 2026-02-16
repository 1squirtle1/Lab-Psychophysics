% Analyze_Adaptation_Staircase.m
clear all;
close all;


% % For checking single data files use this format - easeir than importing
% % .mat file list
% datalist = ["102423_oo_6_Data.mat"];

% Define the file containing the list of .mat file names
cd P:\labFolderNew\Steven\MotionRepulsionData\Adaptation\PsychophysicsData\Owen_data\matFiles
[fileList,pname] = uigetfile('*.txt','pick a filename file');

% Read the file names from the text file
fileID = fopen(fileList, 'r');
fileNames = textscan(fileID, '%s', 'Delimiter', '\n');
fileNames = fileNames{1};
fclose(fileID);

cd P:\labFolderNew\Steven\MotionRepulsionData\Adaptation\PsychophysicsData\Owen_data

% Initialize a cell array to store the loaded data
data = cell(1, numel(fileNames));

% Loop through the file names and load the corresponding .mat files
for i = 1:numel(fileNames)
    try
        data{i} = load(fileNames{i}); % Load the .mat file
    catch
        fprintf('Error loading %s\n', fileNames{i});
    end

% Access the loaded data, e.g., data{1}, data{2}, etc.
trialdata = data{i}(1);

responsedata = cell2mat({trialdata.adaptation_data.response}); 
degreesdata = cell2mat({trialdata.adaptation_data.degrees}); 

% Initialize variables
    reversals = [];
    valuesAtReversals = [];
    averageReversals = [];
% Define the threshold for a reversal
    threshold = 2;

% % Loop through each data cell
% for i = 1:numel(responsedata)

    % Initialize variables for this data set
    upCount = 0;
    downCount = 0;
    lastTwoResponses = zeros(1, 2);

    % Loop through responses in this data set
    for j = 1:numel(responsedata)
        if responsedata(j) == 1
            upCount = upCount + 1;
        else
            downCount = downCount + 1;
        end
     % Store the last two responses     
        lastTwoResponses(mod(j, 2) + 1) = responsedata(j);
        
     % Check for reversals
        if upCount >= threshold && downCount >= threshold
        if lastTwoResponses(2) == lastTwoResponses(1)
            reversals = [reversals, j];
            if responsedata(j) == 1
               downCount = 0;
            else
                upCount = 0;
            end
        end
        end
     
        
    end
   
% Find the degree values at each reversal
    valuesAtReversals = degreesdata(sub2ind(size(degreesdata), reversals));

% Calculate the average of the last 3 responses
    averageReversals = mean(valuesAtReversals(:,5:end));

% Plot the staircase for this data set
    figure;
    plot(1:numel(responsedata), degreesdata, '-o');
    hold on;
    plot(reversals,valuesAtReversals,'r*')
    title(['Staircase for ', fileNames{i},':',' Avg = ',num2str(averageReversals)]);
    xlabel('Trial');
    ylabel('Degree Separation');
 
% Display reversals and average of the last 3
    disp('Reversals:');
    disp(valuesAtReversals);
    disp('Average of the Last 3:');
    disp(averageReversals);

reversals_bytrial(:,i) = valuesAtReversals;
avgreversals_bytrial(:,i) = averageReversals;
    
end    