classdef PowerGrid < handle
    %PowerGrid class capturing electrical power grid.
    %   Detailed explanation goes here

    properties
        casename
        mpc
        baseVal
        sensorTable
        TrueGridState
        TruePowerDem
        TruePowerGen
        TruePowerFlow
        TrueNetPower
        Bus
        Ybus
        Yf
        Yt
        Bp
        Bpp

    end

    methods
        %% Constructor
        function obj = PowerGrid(casename)
            %PowerGrid Construct an instance of this class
            %   Detailed explanation goes here

            obj.casename = casename;
            obj.mpc = obj.loadMatPowerCase();
            obj.baseVal = obj.mpc.baseMVA;

            [obj.Ybus, obj.Yf, obj.Yt, obj.Bp, obj.Bpp] = getGridMat(obj);

            PowerFlowResults = obj.runPowerFlow();

            obj.Bus = obj.collectBusData(PowerFlowResults);
        end

        %% Getters
        function nBus = getNumBus(obj)
                
            nBus = numel(obj.Bus);
        end

        function baseVal = getBaseVal(obj)
            
            baseVal = obj.baseVal;
        end
        function trueval = getTrueVal(obj, iBus, reqDataType)

            % Get true values table
            BusTrueVals = obj.Bus(iBus).trueVal;

            % Get available measurements data types
            avalDataTypes = string(BusTrueVals.Properties.VariableNames);

            if all(ismember(reqDataType, avalDataTypes))
                trueval =  BusTrueVals{:,reqDataType};
            else
                unAvalDataTypeIdx = ~ismember(reqDataType, avalDataTypes);
                unAvalDataType = reqDataType(unAvalDataTypeIdx);

                fprintf('Measurement data type %s is/are not available at bus %d\n', char(strjoin(unAvalDataType, ', ')), iBus);
                fprintf('Available data is\n');
                disp(BusTrueVals);
                error('Measurement data type mismatch!')
            end


        end

        
        function fullScale = getFullScale(obj, iBus, reqDataType)
            
            BusMeasFullScale = obj.Bus(iBus).fullscale;

            % Get available measurements data types
            avalDataTypes = string(BusMeasFullScale.Properties.VariableNames);

            % Check if requested data type is available
            if all(ismember(reqDataType, avalDataTypes))
                fullScale =  BusMeasFullScale{:,reqDataType};
            else
                unAvalDataTypeIdx = ~ismember(reqDataType, avalDataTypes);
                unAvalDataType = reqDataType(unAvalDataTypeIdx);

                fprintf('Measurement full-scale data type %s is/are not available at bus %d\n', char(strjoin(unAvalDataType, ', ')), iBus);
                fprintf('Available data is\n');
                disp(BusTrueVals);
                error('Measurement full-scale data type mismatch!')
            end

        end

         function [Ybus, Yf, Yt, Bp, Bpp] = getGridMat(obj)

            [Ybus, Yf, Yt] = makeYbus(obj.mpc);
            [Bp, Bpp] = makeB(obj.mpc, 'FDXB');

        end

        function busType = getBusType(obj, busIdx)

            define_constants;
            % 1. Get the bus type column index from MATPOWER constants
            % This ensures compatibility even if MATPOWER updates its structure
            [~, BUS_TYPE] = idx_bus;

            % 2. Extract the raw numeric types (1, 2, 3, or 4)
            rawType = obj.mpc.bus(busIdx, BUS_TYPE);

            % 3. Define the mapping
            % Note: MATPOWER uses 1=PQ, 2=PV, 3=Ref, 4=Isolated
            map = {'PQ', 'PV', 'Ref', 'None'};

            % 4. Get bus type
            busType = map(rawType);
        end

        function measSelector = getMeasSel(obj)
            
            measSelector = obj.sensorTable;
        end
        function Bp = getBp(obj)
            
            Bp = obj.Bp;
        end

        function Bpp = getBpp(obj)
            
            Bpp = obj.Bpp;
        end
        %%
        function mpc = loadMatPowerCase(obj)
            %loadMatPowerCase Summary of this method goes here
            %   Detailed explanation goes here

            mpc = loadcase(obj.casename);

        end

        function PowerFlowResults = runPowerFlow(obj)

            opt = mpoption('verbose', 0, 'out.all', 0);
            define_constants;
            PowerFlowResults = runpf(obj.mpc, opt);

            % Get net power from MATPOWER function
            S_inj_true = makeSbus(obj.baseVal, PowerFlowResults.bus, ...
                PowerFlowResults.gen) * obj.baseVal;
            % S_inj_true = makeSbus(obj.baseVal, obj.mpc.bus,...
            %                       obj.mpc.gen) * obj.baseVal;

            % Store true states (voltage mag. and angle)
            obj.TrueGridState = table( PowerFlowResults.bus(:,BUS_I),  ...
                PowerFlowResults.bus(:,VM),  ...
                PowerFlowResults.bus(:,VA), ...
                'VariableNames', {'Bus','VM','VA'});

            % Store true power generation
            % !! Look into how to compile a table
            obj.TruePowerGen = table( PowerFlowResults.gen(:, GEN_BUS), ...
                PowerFlowResults.gen(:, PG),  ...
                PowerFlowResults.gen(:, QG), ...
                'VariableNames', {'Bus','PG','QG'});

            % Store true power injection
            obj.TruePowerDem = table( PowerFlowResults.bus(:,BUS_I), ...
                PowerFlowResults.bus(:,PD),  ...
                PowerFlowResults.bus(:,QD), ...
                'VariableNames', {'Bus','PD','QD'});

            % Store true power flow
            obj.TruePowerFlow = table(PowerFlowResults.branch(:,1),  ...
                PowerFlowResults.branch(:,2),  ...
                PowerFlowResults.branch(:,PF), ...
                PowerFlowResults.branch(:,QF), ...
                PowerFlowResults.branch(:,PT), ...
                PowerFlowResults.branch(:,QT), ...
                'VariableNames', {'From','To','PF','QF','PT','QT'});


            % !!
            obj.TrueNetPower = table(PowerFlowResults.bus(:,BUS_I), ...
                real(S_inj_true), ... % net PD
                imag(S_inj_true), ... % net QD
                'VariableNames', {'Bus','PD','QD'});


        end
       

        function Bus = collectBusData(obj, PowerFlowResults)
            define_constants;
            nBus = size(obj.mpc.bus, 1);

            % Collect VM, PD, QD
            busMat = obj.mpc.bus(:, [VM, PD, QD]);

            % Extract Full-Scale Values for VM, PD, QD
            busMat_fs = [obj.mpc.bus(:, VMAX), abs(obj.mpc.bus(:, PD)), ...
                abs(obj.mpc.bus(:, QD))];

            % Collect PG, QG from the SOLVED results
            genMat = zeros(nBus, 2);
            genMat_fs = zeros(nBus, 2);

            for iBus = 1:nBus
                % Use PowerFlowResults.gen to get the post-simulation true states
                if ismember(iBus, PowerFlowResults.gen(:, GEN_BUS))
                    idx = (iBus == PowerFlowResults.gen(:, GEN_BUS));
                    genMat(iBus,:) = sum(PowerFlowResults.gen(idx, [PG, QG]), 1);

                    % Full-Scale Rating based on machine limits
                    p_fs = max(abs(obj.mpc.gen(idx, PMAX)), abs(obj.mpc.gen(idx, PMIN)));
                    q_fs = max(abs(obj.mpc.gen(idx, QMAX)), abs(obj.mpc.gen(idx, QMIN)));
                    genMat_fs(iBus,:) = [sum(p_fs), sum(q_fs)];

                else
                    genMat(iBus,:) = [0, 0];
                    genMat_fs(iBus,:) = [0 0];
                end
            end

            % Populate bus property
            for iBus = 1:nBus
                Bus(iBus).id = iBus;
                Bus(iBus).Type = obj.getBusType(iBus);

                % --- FIXED: STRICTLY F_BUS ONLY ---
                % Only find branches where this bus is the Sending End
                branchIdx = find(obj.mpc.branch(:, F_BUS) == iBus);

                % --- EXTRACT FLOW AND CREATE HEADERS DYNAMICALLY ---
                flowData = [];
                flowData_fs = [];
                flowNames = {};

                if ~isempty(branchIdx)
                    % Loop through all valid outgoing branches connected to this bus
                    for k = 1:length(branchIdx)
                        br = branchIdx(k);

                        % Extract strictly PF and QF
                        bus_PF = PowerFlowResults.branch(br, PF);
                        bus_QF = PowerFlowResults.branch(br, QF);

                        % Append to the arrays
                        flowData = [flowData, bus_PF, bus_QF];

                        % Extract Full-Scale Line Rating
                        rateA = obj.mpc.branch(br, RATE_A);
                        if rateA == 0
                            rateA = obj.mpc.baseMVA; 
                        end
                        flowData_fs = [flowData_fs, rateA, rateA];

                        % Create the headers using the actual branch ID
                        flowNames = [flowNames, {char("PF-" + string(br)), char("QF-" + string(br))}];
                    end
                else
                    % Leave data empty
                    
                end

                % Assemble the true value table dynamically
                % Order: VM, PD, QD, [All PF/QF Flows], PG, QG
                Bus(iBus).trueVal = array2table([busMat(iBus, :), flowData, genMat(iBus, :)], ...
                    'VariableNames', [{'VM', 'PD', 'QD'}, flowNames, {'PG', 'QG'}]);

                
                % Assemble the Full-Scale Table matching the exact structure
                Bus(iBus).fullscale = array2table([busMat_fs(iBus, :), flowData_fs, genMat_fs(iBus, :)], ...
                    'VariableNames', [{'VM', 'PD', 'QD'}, flowNames, {'PG', 'QG'}]);
             
                Bus(iBus).meas = [];
            end
        end

        function Bus = collectBusDataWorking(obj, PowerFlowResults)

            define_constants;
            % mpc = obj.loadMatPowerCase(obj.casename);
            % mpc = obj.mpc;
            % PowerFlowResults = obj.runPowerFlow();

            nBus = obj.mpc.bus(end,1);
            % Collect VM, PD, QD
            busMat = obj.mpc.bus(:, [VM, PD, QD]);
            % Collect PF, QF
            % branchMat = PowerFlowResults.branch(:, [PF QF]);
            % branchMat = PowerFlowResults.branch(:, [T_BUS PF QF]);

            % Collect PG, QG
            genMat = zeros(nBus, 2);
            for iBus = 1:nBus
                if ismember(iBus, obj.mpc.gen(:, GEN_BUS))
                    idx = (iBus == obj.mpc.gen(:, GEN_BUS));
                    % In the case of having multiple generators
                    genMat(iBus,:) = sum(obj.mpc.gen(idx,[PG QG]),1);
                    % genMat(iBus,:) = obj.mpc.gen(idx,[PG QG]);
                else
                    genMat(iBus,:) = [0 0];
                end
            end

            % Populate bus property
            for iBus = 1:nBus
                Bus(iBus).id = iBus;
                Bus(iBus).Type = obj.getBusType(iBus);

                % --- FIND THE CORRECT BRANCH FOR THIS BUS ---
                % Look for an outgoing branch first
                % branchIdx = find(obj.mpc.branch(:, F_BUS) == iBus, 1);
                % Include all branches
                branchIdx = find(obj.mpc.branch(:, F_BUS) == iBus);

                % If no outgoing branch, look for an incoming branch
                if isempty(branchIdx)
                    branchIdx = find(obj.mpc.branch(:, T_BUS) == iBus, 1);
                end

                % --- EXTRACT FLOW AND CREATE HEADERS DYNAMICALLY ---
                flowData = [];
                flowNames = {};

                if ~isempty(branchIdx)
                    % Loop through all branches connected to this bus
                    for k = 1:length(branchIdx)
                        br = branchIdx(k);

                        % Get the true power flow values for this specific branch
                        bus_PF = PowerFlowResults.branch(br, PF);
                        bus_QF = PowerFlowResults.branch(br, QF);

                        % Append to the arrays
                        flowData = [flowData, bus_PF, bus_QF];

                        % Create the headers using the actual branch ID (row index)
                        flowNames = [flowNames, {char("PF-" + string(br)), char("QF-" + string(br))}];
                    end
                else
                    % Fallback for an isolated bus (rare, but prevents crashes)
                    flowData = [0, 0];
                    flowNames = {'PF-NaN', 'QF-NaN'};
                end

                % Assemble the true value table dynamically
                % Order: VM, PD, QD, [All PF/QF Flows], PG, QG
                Bus(iBus).trueVal = array2table([busMat(iBus, :), flowData, genMat(iBus, :)], ...
                    'VariableNames', [{'VM', 'PD', 'QD'}, flowNames, {'PG', 'QG'}]);
   
                Bus(iBus).meas = [];
            end

        end

      function createMeasStrategy(obj, caseName)
            % CREATEMEASSTRATEGY Generates a measurement plan based on the specific case.
            
            define_constants;
            nb = size(obj.mpc.bus, 1);
            
            % Initialize lists
            busList = [];     
            branchList = [];  
            % Initialize typeList as a string array to prevent data type mismatch crashes
            typeList = strings(0, 1);    
            usedBranches = [];
            
            % 1. Pre-calculate Total Generation per Bus
            if isfield(obj.mpc, 'gen') && ~isempty(obj.mpc.gen)
                P_Gen_Total = sparse(obj.mpc.gen(:, GEN_BUS), 1, obj.mpc.gen(:, PG), nb, 1);
            else
                P_Gen_Total = zeros(nb, 1);
            end

            % Ensure caseName is a lowercase string for reliable matching
            switch lower(string(caseName))
                % =========================================================
                % IEEE 14-BUS LOGIC (FIXED FOR F_BUS OBSERVABILITY)
                % =========================================================
                case {"case14", "ieee14", "case14.m"}
                    
                    for i = 1:nb
                        % Assign VM to Generator Buses
                        if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                            busList = [busList; i];
                            branchList = [branchList; NaN];
                            typeList = [typeList; "VM"];
                        end
                    end
                    
                    for i = 1:nb
                        % Injections
                        if obj.mpc.bus(i, BUS_TYPE) >= 2 
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PG"; "QG"];
                        end
                        % Loads
                        if abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PD"; "QD"];
                        end
                    end
                    
                    % --- CORRECTED FLOW ASSIGNMENT ---
                    % We only assign flow measurements to the SENDING END (F_BUS) 
                    % to match the collectBusData extraction logic.
                    for i = 1:nb
                        % Find branches where this bus is the sender
                        outBranches = find(obj.mpc.branch(:, F_BUS) == i);
                        
                        % Bus 7 is a major junction in the 14-bus case, 
                        % we allow it to measure more than one branch.
                        if i == 7
                            targetMeas = 2; 
                        else
                            targetMeas = 1; 
                        end
                        
                        measCount = 0;
                        for k = 1:length(outBranches)
                            brIdx = outBranches(k);
                            
                            % Only add if the branch hasn't been instrumented yet
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

                % =========================================================
                % IEEE 30-BUS LOGIC (ENHANCED HUB STRATEGY)
                % =========================================================
                case {"case30", "ieee30", "case30.m", "case_ieee30", "case_ieee30.m"}
                    
                    % Define primary RTU Hubs (For flow measurements)
                    hub_buses = [1, 2, 5, 6, 8, 9, 10, 12, 15, 22, 25, 27];
                    
                    % 1. VM Sensors (Ubiquitous monitoring)
                    for i = 1:nb
                        busList = [busList; i];
                        branchList = [branchList; NaN];
                        typeList = [typeList; "VM"];
                    end
                    
                    % 2. Injection Sensors (PG/QG, PD/QD)
                    for i = 1:nb
                        % PG/QG Sensors: Assigned to Generators/Condensers, excluding 5, 8, and 13.
                        % Bus 11 is physically instrumented here.
                        if (obj.mpc.bus(i, BUS_TYPE) >= 2 || ismember(i, [11, 13])) && ~ismember(i, [5, 8, 13])
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PG"; "QG"];
                        end
                        
                        % PD/QD Pseudo-measurements:
                        % EXCLUSION: Buses 11 and 28 are removed to prevent GMM training failures.
                        % Bus 5 remains to ensure its zero-injection baseline is captured.
                        if (abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4 || i == 5) && ...
                           ~ismember(i, [11, 28])
                            busList = [busList; i; i];
                            branchList = [branchList; NaN; NaN];
                            typeList = [typeList; "PD"; "QD"];
                        end
                    end
                    
                    % 3. Flow Sensors (PF, QF) - Strictly F_BUS only per your Hub definition
                    for i = 1:length(hub_buses)
                        hub = hub_buses(i);
                        candidates = find(obj.mpc.branch(:, F_BUS) == hub);
                        
                        for k = 1:length(candidates)
                            brIdx = candidates(k);
                            
                            % EXCEPTION: Remove Branch 13 flows to simplify Bus 9
                            if brIdx == 13
                                continue; 
                            end
                            
                            busList = [busList; hub; hub];
                            branchList = [branchList; brIdx; brIdx];
                            typeList = [typeList; "PF"; "QF"];
                        end
                    end
                    
                    % 4. Manual Overrides
                    % Explicitly add Branch 36 flows to Bus 28
                    busList = [busList; 28; 28];
                    branchList = [branchList; 36; 36];
                    typeList = [typeList; "PF"; "QF"];
                    
                otherwise
                    error('Measurement strategy not explicitly defined for case: %s', caseName);
            end

            %% Create and Sort Output Table
            sensorTable = table(busList, typeList, branchList, ...
                'VariableNames', {'Bus', 'DataType', 'Branch'});
            
            % Sort by Bus ID for a clean logical summary
            if height(sensorTable) > 0
                sensorTable = sortrows(sensorTable, 'Bus');
            end
            
            % Store the generated strategy in the class property
            obj.sensorTable = sensorTable;
      end
      
      function createMeasStrategyWorking16bus(obj)
            % CREATEMEASSTRATEGY Generates a measurement plan.
            % UPDATED: Removed 'unique' to preserve contiguous PG/QG and PF/QF pairing.

            define_constants;
            nb = size(obj.mpc.bus, 1);

            % Initialize lists
            busList = [];     % Location of the sensor
            branchList = [];  % ID of the branch (NaN for Node sensors)
            typeList = [];    % Measurement Type

            % Track usage to avoid duplicate measurements on the same branch
            usedBranches = [];

            % 1. Pre-calculate Total Generation per Bus
            if isfield(obj.mpc, 'gen') && ~isempty(obj.mpc.gen)
                P_Gen_Total = sparse(obj.mpc.gen(:, GEN_BUS), 1, obj.mpc.gen(:, PG), nb, 1);
            else
                P_Gen_Total = zeros(nb, 1);
            end

            %% 1. Voltage Sensors (VM)
            for i = 1:nb
                % Measure only if Slack (3) or PV (2).
                if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                    busList = [busList; i];
                    branchList = [branchList; NaN];
                    typeList = [typeList; "VM"];
                end
            end

            %% 2. Injection Sensors (PG/QG, PD/QD)
            for i = 1:nb
                % --- A. Generators ---
                if obj.mpc.bus(i, BUS_TYPE) >= 2 && abs(P_Gen_Total(i)) > 1e-3
                    busList = [busList; i; i];
                    branchList = [branchList; NaN; NaN];
                    typeList = [typeList; "PG"; "QG"];
                end

                % --- B. Loads ---
                if abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4
                    busList = [busList; i; i];
                    branchList = [branchList; NaN; NaN];
                    typeList = [typeList; "PD"; "QD"];
                end
            end

            %% 3. Flow Sensors (PF, QF)
            for i = 1:nb
                % Determine how many flow measurements we want for this bus
                if i == 7
                    targetMeas = 2; % Redundancy for Bus 7
                else
                    targetMeas = 1; % Standard
                end

                measCount = 0;

                % Gather all candidate branches connected to this bus
                outBranches = find(obj.mpc.branch(:, F_BUS) == i);
                inBranches  = find(obj.mpc.branch(:, T_BUS) == i);
                candidates = [outBranches; inBranches];

                % Iterate through candidates until target is met
                for k = 1:length(candidates)
                    brIdx = candidates(k);

                    % FORCED: Bypass the usedBranches check ONLY for Bus 8
                    if ~ismember(brIdx, usedBranches) || i == 8
                        % Add Measurement
                        busList = [busList; i; i];
                        branchList = [branchList; brIdx; brIdx];
                        typeList = [typeList; "PF"; "QF"];

                        % Mark as used
                        usedBranches = [usedBranches; brIdx];
                        measCount = measCount + 1;
                    end

                    if measCount >= targetMeas
                        break;
                    end
                end
            end

            %% 4. Create and Sort Output Table
            sensorTable = table(busList, typeList, branchList, ...
                'VariableNames', {'Bus', 'DataType', 'Branch'});

            % Sort strictly by Bus ID. This is a stable sort that preserves
            % the logical block order of the measurements.
            sensorTable = sortrows(sensorTable, 'Bus');

            % Store in class
            obj.sensorTable = sensorTable;
        end


        %% Helper methods
        function printMeasSummary(obj)
            % PRINTMEASSUMMARY Prints a formatted table of all measurements
            % collected at each bus, including branch IDs for flow sensors.
            % Prints all buses in the grid, flagging unmonitored ones.

            if isempty(obj.sensorTable)
                fprintf('Sensor table is empty. Please run createMeasStrategy first.\n');
                return;
            end

            casename = obj.casename;
            tbl = obj.sensorTable; % Measurements table
            nb = size(obj.mpc.bus, 1); % Total number of physical buses

            fprintf('\n=========================================================\n');
            fprintf(['                       ', casename,'                    \n']); 
            fprintf('               MEASUREMENT STRATEGY SUMMARY              \n');
            fprintf('=========================================================\n');
            fprintf(' Bus | Collected Measurements\n');
            fprintf('---------------------------------------------------------\n');

            for b = 1:nb
                % Find all rows in the table for this specific bus
                idx = find(tbl.Bus == b);

                if isempty(idx)
                    % Flag buses that have absolutely zero measurements
                    fprintf(' %3d | *** [NONE - UNMONITORED] ***\n', b);
                else
                    measStrs = strings(length(idx), 1);

                    for k = 1:length(idx)
                        mType = string(tbl.DataType(idx(k)));
                        brID = tbl.Branch(idx(k));

                        % If it is a flow measurement, append the branch ID
                        if (mType == "PF" || mType == "QF") && ~isnan(brID)
                            measStrs(k) = sprintf('%s(br%d)', mType, brID);
                        else
                            measStrs(k) = mType;
                        end
                    end

                    % Join the strings with a comma for clean reading
                    measList = strjoin(measStrs, ', ');

                    % Print the row
                    fprintf(' %3d | %s\n', b, measList);
                end
            end

            fprintf('---------------------------------------------------------\n');
            fprintf(' Total Measurements: %d\n', height(tbl));
            fprintf('=========================================================\n\n');
        end

        function busSensorSelector = table2SensorStruct(obj, sensorTable)
            % TABLE2SENSORSTRUCT Converts the modern sensorTable to the legacy struct array format.
            %
            % Inputs:
            %   sensorTable : Table with columns {ID, Type}
            %
            % Output:
            %   busSensorSelector : nb x 1 struct array with fields:
            %                       (VM, PD, QD, PF, QF, PG, QG) as booleans.

            % 1. Setup
            define_constants;
            mpc = obj.mpc;
            nb = size(mpc.bus, 1);

            % 2. Initialize the array of structs (All False)
            template = struct('VM',false, 'PD',false, 'QD',false, ...
                'PF',false, 'QF',false, 'PG',false, 'QG',false);
            busSensorSelector = repmat(template, nb, 1);

            % 3. Iterate through the Table
            nMeas = height(sensorTable);

            for k = 1:nMeas
                id = sensorTable.ID(k);
                type = char(sensorTable.Type(k)); % Convert string to char

                switch type
                    case {'VM', 'PD', 'QD', 'PG', 'QG'}
                        % For Nodal measurements, ID corresponds directly to Bus Index
                        busSensorSelector(id).(type) = true;

                    case {'PF', 'QF'}
                        % For Flow measurements, ID corresponds to Branch Index.
                        % We must map this Branch ID back to its "From Bus"
                        % to set the flag on the correct bus struct.

                        % Look up From Bus in the branch matrix
                        fromBus = mpc.branch(id, F_BUS);

                        % Set the flag on the Bus Struct
                        busSensorSelector(fromBus).(type) = true;
                end
            end
        end
    end
end