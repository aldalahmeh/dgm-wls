%% Compute MonteCarloPerformance Class (Batch Processing)
%% Clear
clc; clear all; close all;

%% Global Parameters
% Define all the files you want to process in a cell array
fileNames = {
    'dbscanGmmWLS_p_0_001', ...
    'dbscanGmmWLS_p_0_01', ...
    'dbscanGmmWLS_p_0_025', ...
    'dbscanGmmWLS_p_0_05', ...
    'dbscanGmmWLS_p_0_1'
};

filePath = '.\Results';
figsPath = '.\Figures';

%% Monte Carlo Performance Object
% Initialize the object once outside the loop
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Batch Process Loop
for i = 1:length(fileNames)
    currentFile = fileNames{i};
    fprintf('\n--- Processing File %d of %d: %s ---\n', i, length(fileNames), currentFile);
    
    % 1. Load Data 
    Data = MonteCarloPerformObj.load(currentFile);
    
    % 2. Change Algorithm Names
    Data.Algorithms(2).Name = "DGM-WLS";
    Data.Algorithms(3).Name = "IRWLS";
    
    % 3. Compute Performance & Plot
    % This will utilize the newly updated plotLoadedDataMAE with the 0.001 VM ticks
    MonteCarloPerformObj.plotLoadedDataMAE(Data);
    
    % 4. Save Figures
    % Uncommented and set to use the current file name so they don't overwrite
    MonteCarloPerformObj.saveFigures(figsPath, currentFile, 'eps');
    
    % 5. Clean up open figures
    % CRITICAL: Close all figures after saving so the next iteration's 
    % findall(0, ...) doesn't accidentally grab and format old plots!
    close all; 
end

fprintf('\nBatch processing complete! All %d scenarios plotted and saved to %s.\n', length(fileNames), figsPath);