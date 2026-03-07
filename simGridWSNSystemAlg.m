% simGridWSNSystemAlg.m
%% Simulate Grid, WSN/IoT Edge Node System
% This is a top level file that simulates the power grid, WSN/IoT nodes and
% the edge node.
%% Clean up
clear all; clc; close all;

%% 
dbstop if error
dbstop if warning

%% Global Parameters
fileName = 'testWLS';
resultsPath = '.\Results';

%% DEBUG
% DEBUG_FLAG = true;
DEBUG_FLAG = false;

%% Save Results
SAVE_FLAG = false;

%% Grid Parameters
% casename = 'case14';
casename = 'case_ieee30';

%% Simulation Parameters
if DEBUG_FLAG
    nSim = 4e3;
    percentOfData  = 0.5;
else
    % nSim = 4.1e3;
    nSim = 8e3;
    percentOfData  = 0.5;   
end

nTrainningData = ceil(percentOfData * nSim);
trainningDataRange = 1:nTrainningData;
testDataRange      = nTrainningData+1:nSim;
nTestingData = length(testDataRange);

%% Noise Parameters
% Sensor noise standard deviation (p.u.)
% sensorNoiseStd = table( 0.005, 0.02, 0.03, 0.01, ...
%                      0.02, 0.05, 0.01, ...
%                     'VariableNames',{'VM', 'PD', 'QD', 'PF', 'QF','PG', 'QG'});

sensorNoiseStd = table( 0.005, 0.01, 0.02, 0.01, ...
                        0.02, 0.01, 0.02, ...
                       'VariableNames',{'VM', 'PD', 'QD', 'PF', 'QF','PG', 'QG'});



ImpNoise = struct('prob', 0.1, ... % probability
                  'multiplier', 10 ...  % level 
                  );

% DBSCAN-GMM Parameters
% 'eps', 0.6
dbscanGmmParam = struct('dbscanParams', struct('eps', 0.3, 'MinPts', 5), ...
                        'gmmParams',    struct('NumComponents', 2, ...
                                               'Replicates', 1, ...
                                               'RegularizationValue',1e-9, ...
                                                'CovarianceType', 'diagonal'));

% Communication noise SNR
SNR = 20;

maxWlsIter = 100;
wlsTol = 1e-6;

%% Algorithm Selection Configuration
% Define the algorithms to run. 

% Alhorithm selector
algSelect = struct('Ordinary_WLS', true, ...
                   'MAP_WLS',      true, ...
                   'MMSE_WLS',     false, ...
                   'GM_WLS',       true  ...
                    );


%% Algorithms options
% Ordinary WLS
ordWlsOptions.FilteringOverride = true;
ordWlsOptions.WlsInit = 'flat';
ordWlsOptions.DebugMode = DEBUG_FLAG;
ordWlsOptions.weightMatrix = 'meas-specs';
% ordWlsOptions.weightMatrix = 'unity';
ordWlsOptions.RobustEstimator = 'None';
ordWlsOptions.dataImpute = 'None';

% MAP-WLS (data imputation)
mapFlrdWlsOptions.FilteringOverride = true;
mapFlrdWlsOptions.DebugMode = DEBUG_FLAG;
% mapFlrdWlsOptions.WlsInit = 'flat';
mapFlrdWlsOptions.WlsInit = 'warm';
mapFlrdWlsOptions.weightMatrix = 'meas-est';
mapFlrdWlsOptions.RobustEstimator = 'None';
mapFlrdWlsOptions.dataImpute = 'MAP';

% MMSE-WLS (data imputation)
mmseFlrdWlsOptions = mapFlrdWlsOptions;
mmseFlrdWlsOptions.dataImpute = 'MMSE';

% GM-WLS
k = [1.0, 1.2, 1.345, 1.5, 2.0, 3, 5];
gmWlsOptions = ordWlsOptions;
gmWlsOptions.RobustEstimator = 'Huber';
gmWlsOptions.cHuberVal = 1.345;
% k = 1 X
% k = 1.345 X
% k = 1.5 X
% k = 2 X
% k = 3 X
% k = 5

Algorithms = struct('Name', {}, 'Options', {}, 'LineColor', {});

% 1. Ordinary WLS
if algSelect.Ordinary_WLS
    Algorithms(end+1).Name = "Ordinary WLS";
    Algorithms(end).Options = ordWlsOptions;
    Algorithms(end).LineColor = 'b'; % Blue
end

% 2. MAP-WLS 
if algSelect.MAP_WLS
    Algorithms(end+1).Name = "MAP WLS ";
    Algorithms(end).Options = mapFlrdWlsOptions;
    Algorithms(end).LineColor = '#D95319'; % Orange
end

% 3. MMSE-WLS
if algSelect.MMSE_WLS
    Algorithms(end+1).Name = "MMSE WLS";
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
IEEE14BusGrid.createMeasStrategy(casename);

% Print the sensor measurement startegy
IEEE14BusGrid.printMeasSummary()

%% Create Environment Object
IEEE14BusGridEnv = GridEnvironment(IEEE14BusGrid, sensorNoiseStd, ImpNoise, nSim);

% Generate measurements
IEEE14BusGridEnv.genSensorMeas();

% Get sensor noise standard deviation
gridSensorNoiseStd = IEEE14BusGridEnv.getGridSensorNoiseStd();

%% Create Edge Node Object
gridEdgeNodeObj = GridEdgeNode(IEEE14BusGrid, nTestingData, gridSensorNoiseStd, maxWlsIter, wlsTol);

%% Create Channel Object
% commsChannelObj = CommsChannel(SNR);

%% Create Sensors Objects & Train Model
nBus = numel(IEEE14BusGridEnv.Grid.Bus);

fprintf('----------------------------------\n');
fprintf('\n Start DBSCAN-GMM Training...\n')
for iBus = 1:nBus
    % Create sensor nodes objects
    sensorObj(iBus) = GridBusSensorNode(IEEE14BusGridEnv.Grid, iBus, dbscanGmmParam);

    % Train DBSCAN-GMM models
    sensorObj(iBus).trainChains(trainningDataRange);

    % Visualize chain model
    % sensorObj(iBus).visualizeModels();
end
fprintf('\n Training Complete.\n')
fprintf('------------------------------\n');

%%  Inference
fprintf('----------------------------------\n');
fprintf('\n Start DBSCAN-GMM Filtering...\n')
for iBus = 1:nBus

    % Test model
    sensorObj(iBus).gmFiltering(testDataRange);

    % Transmit filtered data
    txData(iBus) = sensorObj(iBus).transmit();

    % Comms channel    
    % rxData(iBus) = commsChannelObj.addWhiteNoise(txData(iBus));
    % Short-circuit channel
    rxData(iBus) = txData(iBus);

    % Receive data at edge node
    gridEdgeNodeObj.receive(rxData(iBus));
end

fprintf('\n Filtering Complete.\n')
fprintf('------------------------------\n');

%% State Estimation
tStart = tic;
fprintf('\n--- Starting Simulation Batch ---\n');

for i = 1:length(Algorithms)
    fprintf('\nRunning: %s\n', Algorithms(i).Name);
    fprintf('------------------------------\n');
    
    % Execute State Estimation
    [estState, convData] = gridEdgeNodeObj.computeStateEst(Algorithms(i).Options);
    % [estState, convData] = gridEdgeNodeObj.computeStateEstDebug(Algorithms(i).Options);
    
    % Store raw results directly into the struct
    Algorithms(i).StateEst = estState;
    Algorithms(i).ConvData = convData;
end

% Print elapsed time
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

%% 
pool = gcp('nocreate');   % return [] if no pool
if ~isempty(pool)
    delete(pool);         % shuts down all workers in that pool
end


if exist('q','var') && isa(q,'parallel.pool.PollableDataQueue')
    close(q);    % prevents further sends
    delete(q);
end