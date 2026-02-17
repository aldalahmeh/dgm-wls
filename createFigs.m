%% Creat Figures and Save

%% Clear
clc; clear all; close all;

%% Global Parameters
fileName = 'testWLS';
filePath = '.\Results';
figsPath = '.\Figures';

%% Simulation ID
simID = '_Test';

%% Monte Carlo Performance Object
MonteCarloPerformObj = MonteCarloPerformance(filePath);

%% Load Data 
Data = MonteCarloPerformObj.load(fileName);

%% Performance Analysis & Plotting

Algorithms = Data.Algorithms;
TrueGridState = Data.TrueGridState;

% Preallocate lists for plotting functions
maeMetricsList = {};
convMetricsList = {};
legendNames = {};

for i = 1:length(Algorithms)
    % 1. Compute MAE Metrics
    Algorithms(i).MAEMetrics = MonteCarloPerformObj.computeMAE(...
        Algorithms(i).StateEst, TrueGridState);
    
    % 2. Compute Convergence Statistics
    Algorithms(i).ConvMetrics = MonteCarloPerformObj.computeConvStats(...
        Algorithms(i).ConvData);
    
    % 3. Collect for Plotting
    maeMetricsList{end+1} = Algorithms(i).MAEMetrics;
    convMetricsList{end+1} = Algorithms(i).ConvMetrics;
    legendNames{end+1} = Algorithms(i).Name;
end

% --- Dynamic Plotting ---
% We use the spread operator (...) to pass the cells as individual arguments

% Plot MAE Comparison
% Assuming plotCompareMAE takes (metric1, metric2, ..., label1, label2, ...)
% Note: You might need to adjust your plotCompareMAE signature to accept a cell array of metrics
% OR use this syntax to unwrap them:
MonteCarloPerformObj.plotCompareMAE(maeMetricsList{:}, legendNames{:});

% Plot Convergence Statistics
MonteCarloPerformObj.plotConvStats(convMetricsList{:}, legendNames{:});

%% Save the Generated Figures
% You can easily change 'pdf' 'eps' to 'png' or 'jpeg'
fileExtension = 'pdf'; 
MonteCarloPerformObj.saveFigures(figsPath, simID, fileExtension);
