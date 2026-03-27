%% Compute MonteCarloPerformance Class

%% Clear
clc; clear all; close all;

%% Global Parameters

% fileName = 'dbscanGmmWLS_p_0_01';
% fileName = 'dbscanGmmWLS_p_0_025';
% fileName = 'dbscanGmmWLS_p_0_05';
fileName = 'dbscanGmmWLS_p_0_1';

% fileName = 'dbscanGmmWLS_p_0_001';
filePath = '.\Results';
figsPath = '.\Figures';

%% Monte Carlo Performance Object
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Load Data 
Data = MonteCarloPerformObj.load(fileName);

%% Change Algorithm Names
Data.Algorithms(2).Name = "DGM-WLS";
Data.Algorithms(3).Name = "IRWLS";

%% Compute Performance
MonteCarloPerformObj.plotLoadedDataMAE(Data);

%% Save Figures
% MonteCarloPerformObj.saveFigures(figsPath, fileName, 'eps');