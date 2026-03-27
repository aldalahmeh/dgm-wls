classdef PowerGrid < handle
    % POWERGRID - Class representing an electrical power grid.
    %   Constructs a grid object for a given MATPOWER case name. This class 
    %   holds MATPOWER case data, admittance matrices, and results from power 
    %   flow analysis to be used in downstream state estimation algorithms.
    
    properties
        casename         % String name of the MATPOWER case (e.g., 'case30')
        mpc              % MATPOWER case struct containing raw grid data
        baseVal          % Base MVA of the system
        sensorTable      % Table defining the measurement placement strategy
        TrueGridState    % Ground truth voltage magnitudes and angles
        TruePowerDem     % Ground truth active and reactive power demand (loads)
        TruePowerGen     % Ground truth active and reactive power generation
        TruePowerFlow    % Ground truth active and reactive branch power flows
        TrueNetPower     % Ground truth net power injections per bus
        Bus              % Array of structures holding node-level data and true values
        Ybus             % System bus admittance matrix
        Yf               % Branch admittance matrix (from end)
        Yt               % Branch admittance matrix (to end)
        Bp               % Fast Decoupled load flow B' matrix
        Bpp              % Fast Decoupled load flow B'' matrix
    end
    
    methods
        %% Constructor
        function obj = PowerGrid(casename)
            % POWERGRID Constructs a PowerGrid instance, loads the case, and initializes true states.
            obj.casename = casename;
            obj.mpc = obj.loadMatPowerCase();
            obj.baseVal = obj.mpc.baseMVA;
            [obj.Ybus, obj.Yf, obj.Yt, obj.Bp, obj.Bpp] = getGridMat(obj);
            PowerFlowResults = obj.runPowerFlow();
            obj.Bus = obj.collectBusData(PowerFlowResults);
        end
        
        %% Getters
        function nBus = getNumBus(obj)
            % GETNUMBUS Returns the total number of buses in the grid.
            nBus = numel(obj.Bus);
        end
        
        function baseVal = getBaseVal(obj)
            % GETBASEVAL Returns the system Base MVA.
            baseVal = obj.baseVal;
        end
        
        function trueval = getTrueVal(obj, iBus, reqDataType)
            % GETTRUEVAL Retrieves requested ground truth measurement columns for a specific bus.
            BusTrueVals = obj.Bus(iBus).trueVal;
            avalDataTypes = string(BusTrueVals.Properties.VariableNames);
            
            if all(ismember(reqDataType, avalDataTypes))
                trueval = BusTrueVals{:, reqDataType};
            else
                unAvalDataTypeIdx = ~ismember(reqDataType, avalDataTypes);
                unAvalDataType = reqDataType(unAvalDataTypeIdx);
                fprintf('Measurement data type %s is/are not available at bus %d\n', char(strjoin(unAvalDataType, ', ')), iBus);
                fprintf('Available data is:\n');
                disp(BusTrueVals);
                error('Measurement data type mismatch!');
            end
        end
        
        function fullScale = getFullScale(obj, iBus, reqDataType)
            % GETFULLSCALE Retrieves full-scale measurement limits for a specific bus.
            BusMeasFullScale = obj.Bus(iBus).fullscale;
            avalDataTypes = string(BusMeasFullScale.Properties.VariableNames);
            
            if all(ismember(reqDataType, avalDataTypes))
                fullScale = BusMeasFullScale{:, reqDataType};
            else
                unAvalDataTypeIdx = ~ismember(reqDataType, avalDataTypes);
                unAvalDataType = reqDataType(unAvalDataTypeIdx);
                fprintf('Measurement full-scale data type %s is/are not available at bus %d\n', char(strjoin(unAvalDataType, ', ')), iBus);
                fprintf('Available data is:\n');
                disp(BusMeasFullScale);
                error('Measurement full-scale data type mismatch!');
            end
        end
        
        function [Ybus, Yf, Yt, Bp, Bpp] = getGridMat(obj)
            % GETGRIDMAT Generates the admittance matrices for the grid.
            [Ybus, Yf, Yt] = makeYbus(obj.mpc);
            [Bp, Bpp] = makeB(obj.mpc, 'FDXB');
        end
        
        function busType = getBusType(obj, busIdx)
            % GETBUSTYPE Translates the MATPOWER bus type integer to a descriptive string.
            define_constants;
            [~, BUS_TYPE] = idx_bus;
            rawType = obj.mpc.bus(busIdx, BUS_TYPE);
            
            % MATPOWER definition: 1=PQ, 2=PV, 3=Ref, 4=Isolated
            map = {'PQ', 'PV', 'Ref', 'None'};
            busType = map(rawType);
        end
        
        function measSelector = getMeasSel(obj)
            % GETMEASSEL Returns the generated measurement strategy table.
            measSelector = obj.sensorTable;
        end
        
        function Bp = getBp(obj)
            % GETBP Returns the B' matrix for Fast Decoupled state estimation.
            Bp = obj.Bp;
        end
        
        function Bpp = getBpp(obj)
            % GETBPP Returns the B'' matrix for Fast Decoupled state estimation.
            Bpp = obj.Bpp;
        end
        
        function impNoisIndicator = getBusImpNoisInd(obj, iBus)
            % GETBUSIMPNOISIND Returns the impulsive noise indicator flag for a given bus.
            impNoisIndicator = obj.Bus(iBus).impNoisIndicator;
        end
        
        %% Core Methods
        function mpc = loadMatPowerCase(obj)
            % LOADMATPOWERCASE Loads the MATPOWER case file into the workspace.
            mpc = loadcase(obj.casename);
        end
        
        function PowerFlowResults = runPowerFlow(obj)
            % RUNPOWERFLOW Executes a Newton-Raphson AC power flow and extracts true system states.
            opt = mpoption('verbose', 0, 'out.all', 0);
            define_constants;
            PowerFlowResults = runpf(obj.mpc, opt);
            
            % Calculate true net power injections
            S_inj_true = makeSbus(obj.baseVal, PowerFlowResults.bus, PowerFlowResults.gen) * obj.baseVal;
            
            % Store true grid states (Voltage Magnitude and Angle)
            obj.TrueGridState = table(PowerFlowResults.bus(:, BUS_I), PowerFlowResults.bus(:, VM), PowerFlowResults.bus(:, VA), ...
                'VariableNames', {'Bus', 'VM', 'VA'});
            
            % Store true power generation
            obj.TruePowerGen = table(PowerFlowResults.gen(:, GEN_BUS), PowerFlowResults.gen(:, PG), PowerFlowResults.gen(:, QG), ...
                'VariableNames', {'Bus', 'PG', 'QG'});
            
            % Store true power demand (loads)
            obj.TruePowerDem = table(PowerFlowResults.bus(:, BUS_I), PowerFlowResults.bus(:, PD), PowerFlowResults.bus(:, QD), ...
                'VariableNames', {'Bus', 'PD', 'QD'});
            
            % Store true branch power flows
            obj.TruePowerFlow = table(PowerFlowResults.branch(:, 1), PowerFlowResults.branch(:, 2), ...
                PowerFlowResults.branch(:, PF), PowerFlowResults.branch(:, QF), ...
                PowerFlowResults.branch(:, PT), PowerFlowResults.branch(:, QT), ...
                'VariableNames', {'From', 'To', 'PF', 'QF', 'PT', 'QT'});
            
            % Store true net power
            obj.TrueNetPower = table(PowerFlowResults.bus(:, BUS_I), real(S_inj_true), imag(S_inj_true), ...
                'VariableNames', {'Bus', 'PD', 'QD'});
        end
       
        function Bus = collectBusData(obj, PowerFlowResults)
            % COLLECTBUSDATA Aggregates true states and full-scale ratings into a per-bus structure.
            define_constants;
            nBus = size(obj.mpc.bus, 1);
            
            busMat = obj.mpc.bus(:, [VM, PD, QD]);
            busMat_fs = [obj.mpc.bus(:, VMAX), abs(obj.mpc.bus(:, PD)), abs(obj.mpc.bus(:, QD))];
            
            genMat = zeros(nBus, 2);
            genMat_fs = zeros(nBus, 2);
            
            for iBus = 1:nBus
                if ismember(iBus, PowerFlowResults.gen(:, GEN_BUS))
                    idx = (iBus == PowerFlowResults.gen(:, GEN_BUS));
                    genMat(iBus,:) = sum(PowerFlowResults.gen(idx, [PG, QG]), 1);
                    
                    p_fs = max(abs(obj.mpc.gen(idx, PMAX)), abs(obj.mpc.gen(idx, PMIN)));
                    q_fs = max(abs(obj.mpc.gen(idx, QMAX)), abs(obj.mpc.gen(idx, QMIN)));
                    genMat_fs(iBus,:) = [sum(p_fs), sum(q_fs)];
                else
                    genMat(iBus,:) = [0, 0];
                    genMat_fs(iBus,:) = [0 0];
                end
            end
            
            for iBus = 1:nBus
                Bus(iBus).id = iBus;
                Bus(iBus).Type = obj.getBusType(iBus);
                
                % Locate branches where this bus is the sending end (F_BUS)
                branchIdx = find(obj.mpc.branch(:, F_BUS) == iBus);
                
                flowData = [];
                flowData_fs = [];
                flowNames = {};
                
                if ~isempty(branchIdx)
                    for k = 1:length(branchIdx)
                        br = branchIdx(k);
                        bus_PF = PowerFlowResults.branch(br, PF);
                        bus_QF = PowerFlowResults.branch(br, QF);
                        
                        flowData = [flowData, bus_PF, bus_QF];
                        
                        rateA = obj.mpc.branch(br, RATE_A);
                        if rateA == 0
                            rateA = obj.mpc.baseMVA; 
                        end
                        flowData_fs = [flowData_fs, rateA, rateA];
                        
                        flowNames = [flowNames, {char("PF-" + string(br)), char("QF-" + string(br))}];
                    end
                end
                
                % Assemble the true values table for the bus
                Bus(iBus).trueVal = array2table([busMat(iBus, :), flowData, genMat(iBus, :)], ...
                    'VariableNames', [{'VM', 'PD', 'QD'}, flowNames, {'PG', 'QG'}]);
                
                % Assemble the full-scale rating table for the bus
                Bus(iBus).fullscale = array2table([busMat_fs(iBus, :), flowData_fs, genMat_fs(iBus, :)], ...
                    'VariableNames', [{'VM', 'PD', 'QD'}, flowNames, {'PG', 'QG'}]);
             
                Bus(iBus).meas = [];
            end
        end

        function createMeasStrategy(obj, caseName)
            % CREATEMEASSTRATEGY Configures the sensor placement strategy based on the grid topology.
            define_constants;
            nb = size(obj.mpc.bus, 1);
            
            busList = [];     
            branchList = [];  
            typeList = strings(0, 1);    
            usedBranches = [];
            
            switch lower(string(caseName))
                case {"case14", "ieee14", "case14.m"}
                    % --- IEEE 14-BUS STRATEGY ---
                    % Voltage Magnitude (VM) Sensors
                    for i = 1:nb
                        if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                            busList = [busList; i];
                            branchList = [branchList; NaN];
                            typeList = [typeList; "VM"];
                        end
                    end
                    
                    % Injection Sensors (PG/QG, PD/QD)
                    for i = 1:nb
                        if obj.mpc.bus(i, BUS_TYPE) >= 2 
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PG"; "QG"];
                        end
                        if abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PD"; "QD"];
                        end
                    end
                    
                    % Flow Sensors (PF, QF)
                    for i = 1:nb
                        outBranches = find(obj.mpc.branch(:, F_BUS) == i);
                        targetMeas = (i == 7) + 1; % Bus 7 gets 2 measurements, others get 1
                        
                        measCount = 0;
                        for k = 1:length(outBranches)
                            brIdx = outBranches(k);
                            if ~ismember(brIdx, usedBranches)
                                busList = [busList; i; i];
                                branchList = [branchList; brIdx; brIdx];
                                typeList = [typeList; "PF"; "QF"];
                                usedBranches = [usedBranches; brIdx];
                                measCount = measCount + 1;
                            end
                            if measCount >= targetMeas
                                break;
                            end
                        end
                    end

                case {"case30", "ieee30", "case30.m", "case_ieee30", "case_ieee30.m"}
                    % --- IEEE 30-BUS STRATEGY ---
                    hub_buses = [1, 2, 5, 6, 8, 9, 10, 12, 15, 22, 25, 27];
                    
                    % Voltage Magnitude (VM) Sensors
                    for i = 1:nb
                        busList = [busList; i];
                        branchList = [branchList; NaN];
                        typeList = [typeList; "VM"];
                    end
                    
                    % Injection Sensors (PG/QG, PD/QD)
                    for i = 1:nb
                        if (obj.mpc.bus(i, BUS_TYPE) >= 2 || ismember(i, [11, 13])) && ~ismember(i, [5, 8, 13])
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PG"; "QG"];
                        end
                        if (abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4 || i == 5) && ~ismember(i, [11, 28])
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PD"; "QD"];
                        end
                    end
                    
                    % Flow Sensors (PF, QF)
                    for i = 1:length(hub_buses)
                        hub = hub_buses(i);
                        candidates = find(obj.mpc.branch(:, F_BUS) == hub);
                        for k = 1:length(candidates)
                            brIdx = candidates(k);
                            if brIdx == 13
                                continue; % Skip branch 13 to simplify Bus 9
                            end
                            busList = [busList; hub; hub];
                            branchList = [branchList; brIdx; brIdx];
                            typeList = [typeList; "PF"; "QF"];
                        end
                    end
                    
                    % Manual override for Branch 36 on Bus 28
                    busList = [busList; 28; 28];
                    branchList = [branchList; 36; 36];
                    typeList = [typeList; "PF"; "QF"];
                    
                otherwise
                    error('Measurement strategy not explicitly defined for case: %s', caseName);
            end
            
            sensorTable = table(busList, typeList, branchList, 'VariableNames', {'Bus', 'DataType', 'Branch'});
            if height(sensorTable) > 0
                sensorTable = sortrows(sensorTable, 'Bus');
            end
            obj.sensorTable = sensorTable;
        end
        
        %% Helper methods
        function printMeasSummary(obj)
            % PRINTMEASSUMMARY Prints a formatted terminal report of all generated measurements.
            if isempty(obj.sensorTable)
                fprintf('Sensor table is empty. Please run createMeasStrategy first.\n');
                return;
            end
            
            casename = obj.casename;
            tbl = obj.sensorTable;
            nb = size(obj.mpc.bus, 1); 
            
            fprintf('\n=========================================================\n');
            fprintf(['                       ', casename,'                    \n']); 
            fprintf('               MEASUREMENT STRATEGY SUMMARY              \n');
            fprintf('=========================================================\n');
            fprintf(' Bus | Collected Measurements\n');
            fprintf('---------------------------------------------------------\n');
            
            for b = 1:nb
                idx = find(tbl.Bus == b);
                if isempty(idx)
                    fprintf(' %3d | *** [NONE - UNMONITORED] ***\n', b);
                else
                    measStrs = strings(length(idx), 1);
                    for k = 1:length(idx)
                        mType = string(tbl.DataType(idx(k)));
                        brID = tbl.Branch(idx(k));
                        if (mType == "PF" || mType == "QF") && ~isnan(brID)
                            measStrs(k) = sprintf('%s(br%d)', mType, brID);
                        else
                            measStrs(k) = mType;
                        end
                    end
                    measList = strjoin(measStrs, ', ');
                    fprintf(' %3d | %s\n', b, measList);
                end
            end
            fprintf('---------------------------------------------------------\n');
            fprintf(' Total Measurements: %d\n', height(tbl));
            fprintf('=========================================================\n\n');
        end
        
        function busSensorSelector = table2SensorStruct(obj, sensorTable)
            % TABLE2SENSORSTRUCT Converts the modern sensor table to a legacy struct array format.
            define_constants;
            mpc = obj.mpc;
            nb = size(mpc.bus, 1);
            
            % Initialize boolean template
            template = struct('VM', false, 'PD', false, 'QD', false, ...
                              'PF', false, 'QF', false, 'PG', false, 'QG', false);
            busSensorSelector = repmat(template, nb, 1);
            nMeas = height(sensorTable);
            
            for k = 1:nMeas
                id = sensorTable.ID(k);
                type = char(sensorTable.Type(k)); 
                
                switch type
                    case {'VM', 'PD', 'QD', 'PG', 'QG'}
                        % Nodal measurements assign directly to Bus Index
                        busSensorSelector(id).(type) = true;
                    case {'PF', 'QF'}
                        % Flow measurements assign to the branch's "From Bus"
                        fromBus = mpc.branch(id, F_BUS);
                        busSensorSelector(fromBus).(type) = true;
                end
            end
        end
    end
end