% simGridWSNSystem.m
%% Simulate Grid, WSN/IoT Edge Node System
% This is a top level file that simulates the power grid, WSN/IoT nodes and
% the edge node.
%% Clean up
clear; clc; close all;

%% 
dbstop if error
dbstop if warning

%% Global Parameters
fileName = 'testWLS';
resultsPath = '.\Results';

%% Simulation Parameters
% nSim = 10e3;
nSim = 4e3;
percentOfData  = 0.2;
nTrainningData = percentOfData * nSim;
trainningDataRange = 1:nTrainningData;
testDataRange      = nTrainningData+1:nSim;
nTestingData = length(testDataRange);

% Sensor noise parameters
sensorNoise = table( 0.005, 0.02, 0.03, 0.01, ...
                     0.02, 0.005, 0.01, ...
                    'VariableNames',{'VM', 'PD', 'QD', 'PF', 'QF','PG', 'QG'});

% Impulsive noise parameters
ImpNoise = struct('prob', 0.05, ... % probability
                  'level', 0.1 ...  % level 
                  );
 
% Communication noise SNR
SNR = 20;

maxWlsIter = 100;
wlsTol = 1e-3;

%% Grid Parameters
casename = 'case14';

%% Create Grid Object
IEEE14BusGrid = PowerGrid(casename);
% Create sensor measurement startegy
IEEE14BusGrid.createMeasStrategyWorking();
% IEEE14BusGrid.createMeasStrategyWeakBus8();

%% Create Environment Object
IEEE14BusGridEnv = GridEnvironment(IEEE14BusGrid, sensorNoise, ImpNoise, nSim);

% Generate measurements
IEEE14BusGridEnv.genSensorMeas();

%% Create Edge Node Object
% gridEdgeNodeObj = GridEdgeNode(IEEE14BusGrid, nTrainningData, maxWlsIter, wlsTol);
gridEdgeNodeObj = GridEdgeNode(IEEE14BusGrid, nTestingData, maxWlsIter, wlsTol);

%% Create Channel Object
commsChannelObj = CommsChannel(SNR);

%% Create Sensors Objects & Train Model
nBus = numel(IEEE14BusGridEnv.Grid.Bus);

for iBus = 1:nBus
    % Create sensor nodes objects
    sensorObj(iBus) = GridBusSensorNode(IEEE14BusGridEnv.Grid, iBus);

    % Train DBSCAN-GMM models
    sensorObj(iBus).trainChains(trainningDataRange);

    % Visualize chain model
    % sensorObj(iBus).visualizeModels();
end


%%  Inference
for iBus = 1:nBus

    % Test model
    sensorObj(iBus).gmFiltering(testDataRange);

    % Transmit filtered data
    txData(iBus) = sensorObj(iBus).transmit();

    % Comms channel
    rxData(iBus) = commsChannelObj.addWhiteNoise(txData(iBus));

    % Receive data at edge node
    gridEdgeNodeObj.receive(rxData(iBus));
end

%% Perform State Estimation
% Ordinary WLS
ordWlsOptions.FilteringOverride = true;
ordWlsOptions.WlsInit = 'flat';
ordWlsOptions.DebugMode = false;
ordWlsOptions.weightMatrix = 'comms';
ordWlsOptions.RobustEstimator = 'None';
ordWlsOptions.dataImpute = 'None';

% Filtered WLS (historical data imputation)
histFlrdWlsOptions.FilteringOverride = true;
histFlrdWlsOptions.WlsInit = 'warm';
histFlrdWlsOptions.DebugMode = false;
histFlrdWlsOptions.weightMatrix = 'meas-comms';
histFlrdWlsOptions.RobustEstimator = 'None';
histFlrdWlsOptions.dataImpute = 'Historical';

% Filtered WLS (MAP data imputation)
mapFlrdWlsOptions = histFlrdWlsOptions;
mapFlrdWlsOptions.dataImpute = 'MAP';

% GM-WLS
gmWlsOptions = ordWlsOptions;
gmWlsOptions.RobustEstimator = 'Huber';

% Compute WLS
tStart = tic;
fprintf('\nOrdinary WLS Simulation\n')
fprintf('------------------------------\n')
[ordWlsStateEst, ordWlsConv] = gridEdgeNodeObj.computeStateEst(ordWlsOptions);

fprintf('\nFiltered WLS (Hist) Simulation\n')
fprintf('------------------------------\n')
[histFltrdWlsStateEst, histFltrdWlsConv] = gridEdgeNodeObj.computeStateEst(histFlrdWlsOptions);

fprintf('\nFiltered WLS (MAP) Simulation\n')
fprintf('------------------------------\n')
[mapFltrdWlsStateEst, mapFltrdWlsConv] = gridEdgeNodeObj.computeStateEst(mapFlrdWlsOptions);

fprintf('\nGM-WLS Simulation\n')
fprintf('------------------------------\n')
[gmdWlsStateEst, gmdWlsConv] = gridEdgeNodeObj.computeStateEst(gmWlsOptions);


totalSec = toc(tStart);
fprintf('computeStateEst total time: %.3f s\n', totalSec);

%% Performance Analysis 
MonteCarloPerformObj = MonteCarloPerformance(resultsPath);

% Compute mean absolute error(MAE)
ordPerfMetrics      = MonteCarloPerformObj.computeMAE(ordWlsStateEst,   IEEE14BusGrid.TrueGridState); 
histFltrPerfMetrics = MonteCarloPerformObj.computeMAE(histFltrdWlsStateEst, IEEE14BusGrid.TrueGridState); 
mapFltrPerfMetrics  = MonteCarloPerformObj.computeMAE(mapFltrdWlsStateEst, IEEE14BusGrid.TrueGridState); 
gmPerfMetrics       = MonteCarloPerformObj.computeMAE(gmdWlsStateEst,   IEEE14BusGrid.TrueGridState); 

% Compute covergence statistics:
%  - Convergence rate.
%  - Average # iteration.

ordWlsConvMetrics      = MonteCarloPerformObj.computeConvStats(ordWlsConv);
histFltrWlsConvMetrics = MonteCarloPerformObj.computeConvStats(histFltrdWlsConv);
mapFltrWlsConvMetrics  = MonteCarloPerformObj.computeConvStats(mapFltrdWlsConv);
gmWlsConvMetrics       = MonteCarloPerformObj.computeConvStats(gmdWlsConv);


%% Plotting
% MonteCarloPerformObj.plotMAE(ordPerfMetrics, 'Ordinary WLS');
% 
% MonteCarloPerformObj.plotMAE(fltrPerfMetrics,'Filtered WLS');

% MonteCarloPerformObj.plotMAE(gmPerfMetrics,'GM-WLS');

% MonteCarloPerformObj.plotCompareMAE(ordPerfMetrics, fltrPerfMetrics, 'Ordinary WLS', 'Filtered WLS')

MonteCarloPerformObj.plotCompareMAE(ordPerfMetrics, histFltrPerfMetrics, mapFltrPerfMetrics, gmPerfMetrics, ...
                                    'Ordinary WLS', 'Filtered WLS (Hist.)', 'Filtered WLS (MAP)', 'GM-WLS')


MonteCarloPerformObj.plotConvStats(ordWlsConvMetrics, histFltrWlsConvMetrics, mapFltrWlsConvMetrics, gmWlsConvMetrics, ...
                                    'Ordinary WLS', 'Filtered WLS (Hist.)', 'Filtered WLS (MAP)', 'GM-WLS')

%% Save results
% MonteCarloPerformObj.save(ordWlsStateEst, fltrdWlsStateEst, IEEE14BusGrid.TrueGridState, fileName);

% MonteCarloPerformObj.computeMAE(wlsStateEst);


