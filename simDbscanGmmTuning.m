% simDbscanGmmTuning.m
%% Tune the DBSCAN-GMM parameters

%% Clean up
clear all; clc; close all;

%% 
dbstop if error
dbstop if warning

%% Global Parameters
fileName = 'tuningDBSCAN_GMM';
resultsPath = '.\Results';

%% DEBUG
% DEBUG_FLAG = true;
DEBUG_FLAG = false;

%% Save Results
% SAVE_FLAG = false;
SAVE_FLAG = true;

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

% DBSCAN-GMM parameters
% 'eps', 0.6
dbscanGmmParam = struct('dbscanParams', struct('eps', 0.3, 'MinPts', 5), ...
                        'gmmParams',    struct('NumComponents', 2, ...
                                               'Replicates', 1, ...
                                               'RegularizationValue',1e-9, ...
                                                'CovarianceType', 'diagonal'));

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

%% Create Sensors Objects Train & Test Model
nBus = numel(IEEE30BusGridEnv.Grid.Bus);

minptsValues  = 3 : 1 : 10;
epsilonValues = 0.1 : 0.05 : 0.4;
epsCnt = 1;
minPtsCnt = 1;

f1Score = zeros(length(minPtsCnt), length(epsilonValues));

fprintf('Performing parameter sweep for the DBSCAN-MinPts & eps\n')

tStart = tic;

for minptsVal = minptsValues
    for epsVal = epsilonValues
        dbscanGmmParam.dbscanParams.eps = epsVal;
        dbscanGmmParam.dbscanParams.MinPts = minptsVal;
        estImpNoisInd = cell(1,nBus);
        fprintf('MinPts = %d  eps = %1.2f\n',minptsVal, epsVal );
        for iBus = 1:nBus
            % Create sensor nodes objects
            sensorObj(iBus) = GridBusSensorNode(IEEE30BusGridEnv.Grid, iBus, dbscanGmmParam);

            % Train DBSCAN-GMM models
            sensorObj(iBus).trainChains(trainningDataRange);

            % Test model
            sensorObj(iBus).gmFiltering(testDataRange);

            % Store corrupt data index over all buses
            estImpNoisInd{iBus} = sensorObj(iBus).getSensorImpNoisInd();

        end        
        f1Score(minPtsCnt, epsCnt) = MonteCarloPerformObj.computeF1Score( trueImpNoisInd, estImpNoisInd);
        fprintf('F1 Score = %1.6f\n', f1Score(minPtsCnt, epsCnt));
        epsCnt = epsCnt + 1;
    end
    epsCnt = 1;
    minPtsCnt = minPtsCnt + 1;
end

% Print elapsed time
totalMin = toc(tStart)/60;
fprintf('\nParameter sweep completed in %.3f min.\n', totalMin);

%% Find best parameters
[M,linIdx] = max(f1Score(:));          % M is max, linIdx is linear index
[row,col] = ind2sub(size(f1Score),linIdx);
% Print the best parameters found
fprintf('\nBest MinPts = %d, Best eps = %1.2f with F1 Score = %1.6f\n', minptsValues(row), epsilonValues(col), M);

%% Save Results
if SAVE_FLAG
% Extract the property into a named workspace variable
    
    % Now pass the clean variable name
    MonteCarloPerformObj.saveF1Score(fileName, minptsValues, epsilonValues, f1Score);
end

%% Create latex table
texFileName = 'F1ScoreTable';
MonteCarloPerformObj.generateF1LatexTable(fileName, texFileName);