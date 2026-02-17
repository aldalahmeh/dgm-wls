%% EdgeNode Test
clc; clear all; close all

%% Instatiate Test Class
t = TestGridEdgeNode();
t.setupEnvironment();

%% Run tests
% t.testNetInjections();
% t.testJacbianMatrix();
% t.testVoltageAngle();
% t.testInitEstWarmZeroInjection();
% t.testInitEstWarmPhysicsP();
t.testWlsConvergenceFlatStart();
