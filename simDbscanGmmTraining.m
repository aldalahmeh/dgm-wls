% simDbscanGmmTraining.m
%% Simulate Grid, WSN/IoT Edge Node System
% This is a top level file that trains the DBSCAN-GMM model at the sensor
% nodes.
%% Clean up
clear all; clc; close all;

%% 
dbstop if error
dbstop if warning

%% Global Parameters
figsPath = '.\Figures';

%% DEBUG
% DEBUG_FLAG = true;
DEBUG_FLAG = false;

%% Save Results
% SAVE_FLAG = false;
SAVE_FLAG = true;

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



ImpNoise = struct('prob', 0.01, ... % probability
                  'multiplier', 10 ...  % level 
                  );

% DBSCAN-GMM Parameters
% 'eps', 0.6
dbscanGmmParam = struct('dbscanParams', struct('eps', 0.2, 'MinPts', 5), ...
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
                   'MME_WLS',      true, ...
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

% MME-WLS (data imputation)
mapFlrdWlsOptions.FilteringOverride = true;
mapFlrdWlsOptions.DebugMode = DEBUG_FLAG;
% mapFlrdWlsOptions.WlsInit = 'flat';
mapFlrdWlsOptions.WlsInit = 'warm';
mapFlrdWlsOptions.weightMatrix = 'meas-est';
mapFlrdWlsOptions.RobustEstimator = 'None';
mapFlrdWlsOptions.dataImpute = 'MME';

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

% 2. MME-WLS 
if algSelect.MME_WLS
    Algorithms(end+1).Name = "MME WLS";
    Algorithms(end).Options = mapFlrdWlsOptions;
    Algorithms(end).LineColor = '#D95319'; % Orange
end

% 3. MMSE-WLS
if algSelect.MMSE_WLS
    Algorithms(end+1).Name = "MMSE WLS";
    Algorithms(end).Options = mmseFlrdWlsOptions;
    Algorithms(end).LineColor = '#EDB120'; % Yellow
end

% 4. GM-WLS (Huber)
if algSelect.GM_WLS
    Algorithms(end+1).Name = "GM-WLS (Huber)";
    Algorithms(end).Options = gmWlsOptions;
    Algorithms(end).LineColor = '#7E2F8E'; % Purple
end

%% Create Grid Object
IEEE30BusGrid = PowerGrid(casename);

% Create sensor measurement startegy
IEEE30BusGrid.createMeasStrategy(casename);

% Print the sensor measurement startegy
IEEE30BusGrid.printMeasSummary()

%% Create Environment Object
IEEE30BusGridEnv = GridEnvironment(IEEE30BusGrid, sensorNoiseStd, ImpNoise, nSim);

% Generate measurements
IEEE30BusGridEnv.genSensorMeas();

% Get sensor noise standard deviation
gridSensorNoiseStd = IEEE30BusGridEnv.getGridSensorNoiseStd();

%% Create Edge Node Object
gridEdgeNodeObj = GridEdgeNode(IEEE30BusGrid, nTestingData, gridSensorNoiseStd, maxWlsIter, wlsTol);

%% Create Channel Object
% commsChannelObj = CommsChannel(SNR);

%% Create Sensors Objects & Train Model
nBus = numel(IEEE30BusGridEnv.Grid.Bus);
candidateBuses = [2, 6, 10, 12, 15];
fprintf('----------------------------------\n');
fprintf('\n Start DBSCAN-GMM Training...\n')
for iBus = 1:nBus
    % Create sensor nodes objects
    sensorObj(iBus) = GridBusSensorNode(IEEE30BusGridEnv.Grid, iBus, dbscanGmmParam);

    % Train DBSCAN-GMM models
    sensorObj(iBus).trainChains(trainningDataRange);

    % Visualize chain model
    if any(iBus == candidateBuses)
        sensorObj(iBus).visualizeModels();
    end
end
fprintf('\n Training Complete.\n')
fprintf('------------------------------\n');


%% Performance Analysis & Plotting
MonteCarloPerformObj = MonteCarloPerformance(resultsPath);

%% 
% Ensure the 2x2 figure is currently active
sourceFig = gcf;

% Define the strict IEEE publication size
targetSize = [100, 100, 400, 320];

% Find all axes, excluding any legends
allAxes = findobj(sourceFig, 'Type', 'axes');
allAxes = allAxes(~ismember(get(allAxes,'Tag'), {'legend'}));

for i = 1:length(allAxes)
    ax = allAxes(i);
    
    % 1. Create a new, strictly sized independent figure
    newFig = figure('Position', targetSize, 'Color', 'w');
    
    % 2. Copy the existing axes into the new figure
    newAx = copyobj(ax, newFig);
    
    % 3. Recreate the legend with exact text and location
    if ~isempty(ax.Legend)
        newLeg = legend(newAx, ax.Legend.String);
        newLeg.Location = ax.Legend.Location;
    end
    
    % 4. THICKEN THE ELLIPSOIDS
    plotLines = findall(newAx, 'Type', 'line');
    for k = 1:length(plotLines)
        if ~strcmpi(plotLines(k).LineStyle, 'none')
            plotLines(k).LineWidth = 2.5; 
        end
    end
    
    % 5. Standardize the inner plotting area to prevent margins from shifting
    newAx.Position = [0.15 0.15 0.75 0.75]; 
    
    % 6. Identify the plot type via its title to apply the correct aspect ratio
    titleStr = strjoin(string(newAx.Title.String));
    if contains(titleStr, 'PD-QD') || contains(titleStr, 'PF-QF')
        axis(newAx, 'equal');
    else
        axis(newAx, 'normal');
    end
    
    % 7. Format a clean filename from the title
    safeName = strrep(titleStr, ' ', '_');
    safeName = strrep(safeName, '-', '_');
    if isempty(safeName)
        safeName = sprintf('Plot_%d', i); 
    end
    
    % 8. REMOVE THE TITLE
    title(newAx, ''); 
    
    % 9. Save the perfectly proportioned standalone figure
    MonteCarloPerformObj.saveTargetFigure(newFig, '.\Figures', [char(safeName), '_Bus15.pdf']);
    
    % Close the temporary figure to keep the workspace clean
    close(newFig);
end

fprintf('Successfully split the layout, thickened ellipsoids, removed titles, and saved PDFs.\n');
%% Save Figures
if SAVE_FLAG
    % Target figure number 5 (assuming figure 5 corresponds to Bus 5)
    candBusNum = 15;
    figNum = find(candBusNum == candidateBuses);

    % Specify the exact filename with the .pdf extension
    pdfFileName = sprintf('DBSCAN_GMM_Training_Bus%d.pdf',candBusNum);

    % Call the new method to format and export it
    MonteCarloPerformObj.saveTargetFigure(figNum, '.\Figs', pdfFileName);
end