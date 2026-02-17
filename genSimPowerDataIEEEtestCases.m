% genSimPowerDataIEEEtestCases.m
%% Generate Power Load Simulated Data from IEEE Test Cases
% Generate random demand power values (P&Q) that are uniformly distributed.
% Then find the rest of the grid values by running MATPOWER power flow
% function. Finally, save data to file.
%% Initialization
clear
clc
close all
dbstop if error
dbstop if warning

%%
addpath('C:\Users\sami\Dropbox\Research\Smart Grid\Simulation Files\Library')

%% Flag
SAVE_FLAG = 1;

%% Simulation parameters
nSim = 4e3;

loadMaxVal = 1.1;
loadMinVal = 0.9;

noiseLevel = struct('VmagLevel' , 0.005, ...
    'PdemLevel' , 0.02, ...
    'QdemLevel' , 0.03, ...
    'PflowLevel', 0.01, ...
    'QflowLevel', 0.02, ...
    'PgenLevel' , 0.005, ...
    'QgenLevel' , 0.01, ...
    'impLevel'  , 0.10);

impProb         = 0.01;       % impulsive noise occurance probability

%% Generate True Grid Values

% Load test case
define_constants;
casename = 'case14';
% casename = 'case_ieee30';

mpc = loadcase(casename);
% Quiet options
opt = mpoption('verbose', 0, 'out.all', 0);

% Extract true data from test case
[Pgrid, Qgrid, Vgrid, Ind] = extrctMpcData(casename);

nBus = mpc.bus(end,1);
nLoadBus = length(Ind.PD);

% Generate randomly uniform distributed instances of demand P and Q data
% PlMeas = Pgrid.PD + Pgrid.PD .* (loadMaxVal - loadMinVal) .* rand(nLoadBus,nSim);
% QlMeas = Qgrid.QD + Qgrid.QD .* (loadMaxVal - loadMinVal) .* rand(nLoadBus,nSim);

%% Create Bus Container Structure
measType = {'V','PD','QD','PF','QF','PG','QG'};
BusType = {'Ref', 'PV', 'PQ'};
Bus(nBus) = struct('IEEE_Test_Case', casename, ...
    'id',             [],       ... % Bus ID,
    'Bus_Type',       [],       ...
    'meas_Type',      [],       ...
    'meas', nan*ones(nSim, length(measType)) );



% Bus type for the IEEE 16-bus test case
ieee16BusType = { ...
    'Ref'; ... % Bus 1  (Slack/Reference)
    'PV';  ... % Bus 2  (Generator)
    'PV';  ... % Bus 3  (Synchronous Condenser)
    'PQ';  ... % Bus 4  (Load)
    'PQ';  ... % Bus 5  (Load)
    'PV';  ... % Bus 6  (Synchronous Condenser)
    'PQ';  ... % Bus 7  (Zero Injection / Transit)
    'PV';  ... % Bus 8  (Synchronous Condenser)
    'PQ';  ... % Bus 9  (Load)
    'PQ';  ... % Bus 10 (Load)
    'PQ';  ... % Bus 11 (Load)
    'PQ';  ... % Bus 12 (Load)
    'PQ';  ... % Bus 13 (Load)
    'PQ'   ... % Bus 14 (Load)
    };




noisParam.noiseLevelImp = noiseLevel.impLevel;

% Feed random power demand values
mpc.bus(Ind.PD,PD) = Pgrid.PD;
mpc.bus(Ind.QD,QD) = Qgrid.QD;

% Run power flow
pfResults = runpf(mpc, opt);


disp('Running measured data generation ...')
for iSim = 1:nSim

    % Impulsive noise indicator for the bus
    impNoisBusInd = (rand(nBus,1) < impProb);

    % Impulsive noise indicator for the branch
    for i = 1:length(pfResults.branch(:,1))
        impNoisBranchInd(i) = impNoisBusInd(pfResults.branch(i,1));
    end

    % Collect VM, PD, QD
    busMat = mpc.bus(:, [VM, PD, QD]);

    % Collect PF, QF
    branchMat = pfResults.branch(:, [PF QF]);

    % Collect PG, QG
    for iBus = 1:nBus
        if ismember(iBus, mpc.gen(:, 1))
            idx = (iBus == mpc.gen(:, 1));
            genMat(iBus,:) = mpc.gen(idx,[PG QG]);
        else
            genMat(iBus,:) = [0 0];
        end
    end

    % Generate measurements with BG noise
    busMeas = busMat + busMat .* [noiseLevel.VmagLevel noiseLevel.PdemLevel noiseLevel.QdemLevel].*randn(size(busMat)) ...
        +  impNoisBusInd .* busMat .* noiseLevel.impLevel .* randn(size(busMat));


    branchMeas = branchMat + branchMat .* [noiseLevel.PflowLevel noiseLevel.QflowLevel].*randn(size(branchMat)) ...
        + impNoisBranchInd' .* branchMat .* noiseLevel.impLevel .* randn(size(branchMat));

    genMeas = genMat + genMat .* [noiseLevel.PgenLevel noiseLevel.QgenLevel].*randn(size(genMat)) ...
        +  impNoisBusInd .* genMat .* noiseLevel.impLevel .* randn(size(genMat));

    % Populate bus structure
    for iBus = 1:nBus
        Bus(iBus).IEEE_Test_Case  = casename;
        Bus(iBus).id = iBus;    
    end

    Bus(iBus).Bus_Type =  ieee16BusType{iBus};
    Bus(iBus).meas_Type = measType;
    Bus(iBus).meas(iSim,:) = [busMeas(iBus, :) branchMeas(iBus, :) genMeas(iBus, :)] ;
end

%% Conver Measured Data to Tables
for iBus = 1:nBus
    Bus(iBus).meas = array2table(Bus(iBus).meas, 'VariableNames',measType);
end


%% Save Data to File
if (SAVE_FLAG)
    disp('Saving data to file ...')

    folder = "./Data";
    fileName = ['Power_Data_', casename] ;
    fullFileName = fullfile(folder, fileName);

    save(fullFileName, 'Bus', 'casename', 'nSim', 'mpc');
    disp('Data saved!')

end
