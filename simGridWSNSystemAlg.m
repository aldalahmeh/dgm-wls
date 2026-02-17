% simGridWSNSystemAlg.m
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

%% DEBUG
DEBUG_FLAG = true;

%% Save Results
SAVE_FLAG = true;

%% Simulation Parameters
if DEBUG_FLAG
    nSim = 4e3;
    percentOfData  = 0.5;
else
    nSim = 10e3;
    percentOfData  = 0.2;
end

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

%% Algorithm Selection Configuration
% Define the algorithms to run. 

% Alhorithm selector
algSelect = struct('Ordinary_WLS',      true, ...
                   'Hist_Filtered_WLS', true, ...
                   'MAP_Filtered_WLS',  false, ...
                   'GM_WLS',            false  ...
                    );


%% Algorithms options
% Ordinary WLS
ordWlsOptions.FilteringOverride = true;
ordWlsOptions.WlsInit = 'flat';
ordWlsOptions.DebugMode = DEBUG_FLAG;
ordWlsOptions.weightMatrix = 'comms';
ordWlsOptions.RobustEstimator = 'None';
ordWlsOptions.dataImpute = 'None';

% Filtered WLS (historical data imputation)
histFlrdWlsOptions.FilteringOverride = true;
histFlrdWlsOptions.WlsInit = 'warm';
histFlrdWlsOptions.DebugMode = DEBUG_FLAG;
histFlrdWlsOptions.weightMatrix = 'meas-comms';
histFlrdWlsOptions.RobustEstimator = 'None';
histFlrdWlsOptions.dataImpute = 'Historical';

% Filtered WLS (MAP data imputation)
mapFlrdWlsOptions = histFlrdWlsOptions;
mapFlrdWlsOptions.dataImpute = 'MAP';

% GM-WLS
gmWlsOptions = ordWlsOptions;
gmWlsOptions.RobustEstimator = 'Huber';

Algorithms = struct('Name', {}, 'Options', {}, 'LineColor', {});

% 1. Ordinary WLS
if algSelect.Ordinary_WLS
    Algorithms(end+1).Name = "Ordinary WLS";
    Algorithms(end).Options = ordWlsOptions;
    Algorithms(end).LineColor = 'b'; % Blue
end

% 2. Filtered WLS (Historical / Phantom Anchor)
if algSelect.Hist_Filtered_WLS
    Algorithms(end+1).Name = "Filtered WLS (Hist)";
    Algorithms(end).Options = histFlrdWlsOptions;
    Algorithms(end).LineColor = '#D95319'; % Orange
end

% 3. Filtered WLS (MAP / Hybrid)
if algSelect.MAP_Filtered_WLS
    Algorithms(end+1).Name = "Filtered WLS (MAP)";
    Algorithms(end).Options = mapFlrdWlsOptions;
    Algorithms(end).LineColor = '#EDB120'; % Yellow
end

% 4. GM-WLS (Huber)
if algSelect.GM_WLS
    Algorithms(end+1).Name = "GM-WLS (Huber)";
    Algorithms(end).Options = gmWlsOptions;
    Algorithms(end).LineColor = '#7E2F8E'; % Purple
end

%% Create Grid Object
IEEE14BusGrid = PowerGrid(casename);
% Create sensor measurement startegy
IEEE14BusGrid.createMeasStrategy();

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


%% State Estimation
tStart = tic;
fprintf('\n--- Starting Simulation Batch ---\n');

for i = 1:length(Algorithms)
    fprintf('\nRunning: %s\n', Algorithms(i).Name);
    fprintf('------------------------------\n');
    
    % Execute State Estimation
    [estState, convData] = gridEdgeNodeObj.computeStateEst(Algorithms(i).Options);
    
    % Store raw results directly into the struct
    Algorithms(i).StateEst = estState;
    Algorithms(i).ConvData = convData;
end

% totalSec = toc(tStart);
% fprintf('\nAll simulations completed in %.3f s\n', totalSec);
totalMin = toc(tStart)/60;
fprintf('\nAll simulations completed in %.3f min.\n', totalMin);

%% Performance Analysis & Plotting
MonteCarloPerformObj = MonteCarloPerformance(resultsPath);

% Preallocate lists for plotting functions
maeMetricsList = {};
convMetricsList = {};
legendNames = {};

for i = 1:length(Algorithms)
    % 1. Compute MAE Metrics
    Algorithms(i).MAEMetrics = MonteCarloPerformObj.computeMAE(...
        Algorithms(i).StateEst, IEEE14BusGrid.TrueGridState);
    
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

%% Save Results
if SAVE_FLAG
% Extract the property into a named workspace variable
    TrueGridState = IEEE14BusGrid.TrueGridState;
    
    % Now pass the clean variable name
    MonteCarloPerformObj.save(fileName, algSelect, Algorithms, TrueGridState);
end