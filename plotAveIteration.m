
%% Clear
clc; clear all; close all;
%%
filePath = '.\Results\';
figsPath = '.\Figures';
fileName = 'average_iter_num';

%%
% Define your files and their matching probabilities
matFiles = {'dbscanGmmWLS_p_0_001.mat', ...
            'dbscanGmmWLS_p_0_025.mat', ...
            'dbscanGmmWLS_p_0_05.mat', ...
            'dbscanGmmWLS_p_0_1.mat'};
probs = [0.01, 0.025, 0.05, 0.10];


%% Monte Carlo Performance Object
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Generate the plot
MonteCarloPerformObj.plotAvgIterTrend(matFiles, probs);

%% Save Figures
MonteCarloPerformObj.saveFigures(figsPath, fileName, 'eps');