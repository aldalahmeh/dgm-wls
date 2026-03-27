%% Batch Script: Average MAE Trend Analysis (VM and VA)
clc; clear all; close all;

%% Global Parameters
% Define the files in order of increasing probability
fileNames = {
    'dbscanGmmWLS_p_0_001', ...
    'dbscanGmmWLS_p_0_01', ...
    'dbscanGmmWLS_p_0_025', ...
    'dbscanGmmWLS_p_0_05', ...
    'dbscanGmmWLS_p_0_1'
};

% Define the corresponding probabilities for the x-axis
probabilities = [0.001, 0.01, 0.025, 0.05, 0.1];

filePath = '.\Results';
figsPath = '.\Figures';

%% Initialize Object
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Generate the Trend Plots
fprintf('\n--- Generating Average MAE Trend Plots ---\n');

% Calls the method to process and plot both VM and VA simultaneously
MonteCarloPerformObj.plotAvgMAETrend(fileNames, probabilities);

%% Save the Figures
% The saveFigures method will detect the two open windows, read their 
% titles to distinguish VM from VA, and save them automatically.
MonteCarloPerformObj.saveFigures(figsPath, 'AvgMAETrend', 'eps');

fprintf('\nBatch processing complete!\n');