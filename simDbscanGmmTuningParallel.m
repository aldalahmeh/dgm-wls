% simDbscanGmmTuning.m
%% Tune the DBSCAN-GMM parameters
%% Clean up
clear all; clc; close all;
%% 
dbstop if error
dbstop if warning

%% Initialize Parallel Pool
% Configure for exactly 20 workers
poolobj = gcp('nocreate');
if isempty(poolobj) || poolobj.NumWorkers ~= 20
    delete(gcp('nocreate'));
    parpool(20);
end

%% Global Parameters
fileName = 'tuningDBSCAN_GMM';
resultsPath = '.\Results';
%% DEBUG
% DEBUG_FLAG = true;
DEBUG_FLAG = false;
%% Save Results
SAVE_FLAG = false;
%% Grid Parameters
casename = 'case_ieee30';
%% Simulation Parameters
if DEBUG_FLAG
    nSim = 4e3;
    percentOfData  = 0.5;
else
    nSim = 10e3;
    percentOfData  = 0.5;   
end
nTrainningData = ceil(percentOfData * nSim);
trainningDataRange = 1:nTrainningData;
testDataRange      = nTrainningData+1:nSim;
nTestingData = length(testDataRange);

%% Noise Parameters
sensorNoiseStd = table( 0.005, 0.01, 0.02, 0.01, ...
                        0.02, 0.01, 0.02, ...
                       'VariableNames',{'VM', 'PD', 'QD', 'PF', 'QF','PG', 'QG'});
ImpNoise = struct('prob', 0.1, ... % probability
                  'multiplier', 10 ...  % level 
                  );

%% Create Grid Object
IEEE30BusGrid = PowerGrid(casename);
% Create sensor measurement startegy
IEEE30BusGrid.createMeasStrategy(casename);

%% Create Environment Object
IEEE30BusGridEnv = GridEnvironment(IEEE30BusGrid, sensorNoiseStd, ImpNoise, nSim);
% Generate measurements
IEEE30BusGridEnv.genSensorMeas();
% Return the corrupted measurements in the grid
trueImpNoisInd = IEEE30BusGridEnv.getGridImpNoisInd(testDataRange);

%% Create Performance Analysis Object
MonteCarloPerformObj = MonteCarloPerformance(resultsPath);

%% Setup Parameter Sweep Grid
nBus = numel(IEEE30BusGridEnv.Grid.Bus);
minptsValues  = 3 : 1 : 10;
epsilonValues = 0.1 : 0.05 : 0.4;

% Flatten the 2D sweep into a 1D list of tasks for parfor
[MinPtsGrid, EpsGrid] = ndgrid(minptsValues, epsilonValues);
minPtsList = MinPtsGrid(:);
epsList    = EpsGrid(:);
nTasks     = length(minPtsList);

% Preallocate a 1D array for the results
f1ScoreList = zeros(nTasks, 1);

fprintf('Performing parameter sweep for DBSCAN-MinPts & eps across %d tasks using 20 workers...\n', nTasks);

%% Execute Parallel Parameter Sweep
tStart = tic;

parfor iTask = 1:nTasks
    % Extract current parameters for this specific worker task
    minptsVal = minPtsList(iTask);
    epsVal    = epsList(iTask);
    
    % Build a strictly local parameter struct to avoid broadcast conflicts
    local_dbscanGmmParam = struct('dbscanParams', struct('eps', epsVal, 'MinPts', minptsVal), ...
                                  'gmmParams',    struct('NumComponents', 2, ...
                                                         'Replicates', 1, ...
                                                         'RegularizationValue', 1e-5, ...
                                                         'CovarianceType', 'diagonal'));
    
    % Preallocate local cell array for the estimated indices
    local_estImpNoisInd = cell(1, nBus);
    
    for iBus = 1:nBus
        % Create temporary sensor node object for this bus
        tempSensor = GridBusSensorNode(IEEE30BusGridEnv.Grid, iBus, local_dbscanGmmParam);
        
        % Train and Test models
        tempSensor.trainChains(trainningDataRange);
        tempSensor.gmFiltering(testDataRange);
        
        % Store corrupt data index
        local_estImpNoisInd{iBus} = tempSensor.getSensorImpNoisInd();
    end        
    
    % Compute and store F1 score for this hyperparameter combination
    f1ScoreList(iTask) = MonteCarloPerformObj.computeF1Score(trueImpNoisInd, local_estImpNoisInd);
    
    fprintf('Completed: MinPts = %d, eps = %1.2f | F1 Score = %1.6f\n', minptsVal, epsVal, f1ScoreList(iTask));
end


% Print elapsed time
totalMin = toc(tStart)/60;
fprintf('\nParameter sweep completed in %.3f min.\n', totalMin);

%% Post-Process and Reshape Results
% Reconstruct the 2D matrix matching the original minptsValues and epsilonValues dimensions
f1ScoreMatrix = reshape(f1ScoreList, size(MinPtsGrid));

fprintf('Parameter sweep complete.\n');

%% Find best parameters
[M,linIdx] = max(f1ScoreMatrix(:));          % M is max, linIdx is linear index
[row,col] = ind2sub(size(f1ScoreMatrix),linIdx);
% Print the best parameters found
fprintf('\nBest MinPts = %d, Best eps = %1.2f with F1 Score = %1.6f\n', minptsValues(row), epsilonValues(col), M);

%% Save Results
if SAVE_FLAG
    % Pass the reshaped 2D matrix for structured saving
    MonteCarloPerformObj.save(fileName, minptsValues, epsilonValues, f1ScoreMatrix);
end