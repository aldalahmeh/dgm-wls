%% Test MonteCarloPerformance Class

%% Clear
clc; clear all; close all;

%% Parameters
fileName = 'dummy';
filePath = '.\Results';

%% Create Dummy Data
data = rand(2000, 1);

%% 
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Save 
MonteCarloPerformObj.save(data, fileName);