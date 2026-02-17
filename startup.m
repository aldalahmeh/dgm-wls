% startup.m
% This script runs automatically when you open MATLAB in this project folder.
clc; clear all
disp('Initializing Grid IoT Project Environment...');

% 1. Define the path to your MATPOWER library
% (Based on the path I saw in your genSimPowerDataIEEEtestCases.m file)
% UPDATE THIS PATH if your MATPOWER folder location has changed.
% matpowerPath = 'C:\Users\sami\Dropbox\Research\Smart Grid\Simulation Files\Library';
matpowerPath = '..\..\matpower8.1';

% 2. Check if the folder exists before adding it
if exist(matpowerPath, 'dir')
    % addpath(genpath(...)) adds the folder AND all subfolders
    addpath(genpath(matpowerPath));
    disp(['[OK] Added to path: ' matpowerPath]);
else
    warning(['[FAIL] MATPOWER path not found: ' matpowerPath]);
    disp('Please update the "matpowerPath" variable in startup.m');
end

% 3. Verify that loadcase is now accessible
if exist('loadcase', 'file')
    disp('[SUCCESS] MATPOWER is ready. You can run your simulation.');
else
    disp('[ERROR] "loadcase" function still not found.');
end