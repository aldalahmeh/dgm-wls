%% Compute MonteCarloPerformance Class

%% Clear
clc; clear all; close all;

%% Global Parameters
fileName = 'testWLS';
filePath = '.\Results';

%% Monte Carlo Performance Object
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Load Data 
Data = MonteCarloPerformObj.load(fileName);

%% Compute Performance
performanceMetrics = MonteCarloPerformObj.computeMAE(Data); 
MonteCarloPerformObj.plotMAE(performanceMetrics);