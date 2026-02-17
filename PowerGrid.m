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

        function Bus = collectBusData(obj, PowerFlowResults)

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

                % --- EXTRACT FLOW AND CREATE HEADERS ---
                % if ~isempty(branchIdx)
                %     % Get the true power flow values for this specific branch
                %     bus_PF = PowerFlowResults.branch(branchIdx, PF);
                %     bus_QF = PowerFlowResults.branch(branchIdx, QF);
                %
                %     % Create the headers using the actual branch ID (row index)
                %     PfColHeader = char("PF-" + string(branchIdx));
                %     QfColHeader = char("QF-" + string(branchIdx));
                % else
                %     % Fallback for an isolated bus (rare, but prevents crashes)
                %     bus_PF = 0; bus_QF = 0;
                %     PfColHeader = "PF-NaN"; QfColHeader = "QF-NaN";
                % end
                %
                % % Assemble the true value table
                % Bus(iBus).trueVal = array2table([busMat(iBus, :), bus_PF, bus_QF, genMat(iBus, :)] ...
                %     , 'VariableNames', {'VM', 'PD', 'QD', PfColHeader, QfColHeader, 'PG', 'QG'});
                % Bus(iBus).trueVal = array2table([busMat(iBus, :), branchMat(iBus, :), genMat(iBus, :)] ...
                %     , 'VariableNames', {'VM', 'PD', 'QD', 'PF', 'QF', 'PG', 'QG'});
                Bus(iBus).meas = [];
                % Bus(iBus).meas = array2table([busMat(iBus, :), branchMat(iBus, :), genMat(iBus, :)] ...
                %                  , 'VariableNames', {'VM', 'PD', 'QD', 'PF', 'QF', 'PG', 'QG'});
            end

        end

        function createMeasStrategy(obj)
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
        function createMeasStrategyWeakBus8(obj)
            % CREATEMEASSTRATEGY Generates a measurement plan.
            % UPDATED: Removed VM on Bus 7. Added DOUBLE Flow measurement on Bus 7.

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
                % REMOVED: Explicit Bus 7 check.
                if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                    busList = [busList; i];
                    branchList = [branchList; NaN]; % Not applicable
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
            % Logic: Find branches connected to this bus.
            % UPDATED: Bus 7 attempts to find 2 branches. Others find 1.
            for i = 1:nb

                % Determine how many flow measurements we want for this bus
                if i == 7
                    targetMeas = 2; % Redundancy for Bus 7
                else
                    targetMeas = 1; % Standard
                end

                measCount = 0;

                % Gather all candidate branches connected to this bus
                % Priority 1: Outgoing (From Bus == i)
                outBranches = find(obj.mpc.branch(:, F_BUS) == i);
                % Priority 2: Incoming (To Bus == i)
                inBranches  = find(obj.mpc.branch(:, T_BUS) == i);

                % Combine candidates
                candidates = [outBranches; inBranches];

                % Iterate through candidates until target is met
                for k = 1:length(candidates)
                    brIdx = candidates(k);

                    if ~ismember(brIdx, usedBranches)
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
            % Sort by Bus ID for readability
            sensorTable = sortrows(sensorTable, 'Bus');
            % Store in class
            obj.sensorTable = sensorTable;
        end

        function createMeasStrategyWeakBus7n8(obj)
            % CREATEMEASSTRATEGY Generates a measurement plan.
            % UPDATED: Separates Bus ID (Location) and Branch ID (Element) into two columns.

            define_constants;
            nb = size(obj.mpc.bus, 1);

            % Initialize lists
            busList = [];     % Location of the sensor
            branchList = [];  % ID of the branch (NaN for Node sensors)
            typeList = [];    % Measurement Type

            % Track usage to avoid duplicate measurements on the same branch
            usedBranches = [];

            % 1. Pre-calculate Total Generation per Bus
            P_Gen_Total = sparse(obj.mpc.gen(:, GEN_BUS), 1, obj.mpc.gen(:, PG), nb, 1);

            %% 1. Voltage Sensors (VM)
            for i = 1:nb
                if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                    busList = [busList; i];
                    branchList = [branchList; NaN]; % Not applicable
                    typeList = [typeList; "VM"];
                end
            end

            %% 2. Injection Sensors (PG/QG, PD/QD)
            for i = 1:nb
                % --- A. Generators ---
                if obj.mpc.bus(i, BUS_TYPE) >= 2 && abs(P_Gen_Total(i)) > 1e-3
                    % Add PG
                    busList = [busList; i];
                    branchList = [branchList; NaN];
                    typeList = [typeList; "PG"];

                    % Add QG
                    busList = [busList; i];
                    branchList = [branchList; NaN];
                    typeList = [typeList; "QG"];
                end

                % --- B. Loads ---
                if abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4
                    % Add PD
                    busList = [busList; i];
                    branchList = [branchList; NaN];
                    typeList = [typeList; "PD"];

                    % Add QD
                    busList = [busList; i];
                    branchList = [branchList; NaN];
                    typeList = [typeList; "QD"];
                end
            end

            %% 3. Flow Sensors (PF, QF)
            % Logic: Find ONE branch connected to this bus to measure.
            for i = 1:nb
                % A. Try finding an OUTGOING branch (From Bus == i)
                branchIdx = find(obj.mpc.branch(:, F_BUS) == i, 1);

                if ~isempty(branchIdx) && ~ismember(branchIdx, usedBranches)
                    % Sensor Location: Bus i
                    % Element: branchIdx
                    busList = [busList; i; i];
                    branchList = [branchList; branchIdx; branchIdx];
                    typeList = [typeList; "PF"; "QF"];

                    usedBranches = [usedBranches; branchIdx];
                else
                    % B. If no outgoing, try INCOMING (To Bus == i)
                    branchIdx = find(obj.mpc.branch(:, T_BUS) == i, 1);

                    if ~isempty(branchIdx) && ~ismember(branchIdx, usedBranches)
                        % Sensor Location: Bus i
                        % Element: branchIdx
                        busList = [busList; i; i];
                        branchList = [branchList; branchIdx; branchIdx];
                        typeList = [typeList; "PF"; "QF"];

                        usedBranches = [usedBranches; branchIdx];
                    end
                end

            end

            %% 4. Create and Sort Output Table
            sensorTable = table(busList, typeList, branchList, ...
                'VariableNames', {'Bus', 'DataType', 'Branch'});

            % Sort by Bus ID for readability
            sensorTable = sortrows(sensorTable, 'Bus');
            % Store in class
            obj.sensorTable = sensorTable;
        end

        function sensorTable = createMeasStrategyWrong(obj)
            % CREATEMEASSTRATEGY Generates a measurement plan.
            % UPDATED: Uses 'unique' to remove duplicate Branch sensors automatically.

            define_constants;
            nb = size(obj.mpc.bus, 1);

            idList = [];
            typeList = [];

            % 1. Pre-calculate Total Generation per Bus
            P_Gen_Total = sparse(obj.mpc.gen(:, GEN_BUS), 1, obj.mpc.gen(:, PG), nb, 1);

            %% 1. Voltage Sensors (VM)
            for i = 1:nb
                % Measure VM at all Voltage-Controlled buses (Slack & PV)
                if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                    idList = [idList; i];
                    typeList = [typeList; "VM"];
                end
            end

            %% 2. Injection Sensors (PG/QG, PD/QD)
            for i = 1:nb
                % --- A. Generators ---
                % Filter out Sync Condensers (PG ~ 0)
                if obj.mpc.bus(i, BUS_TYPE) >= 2 && abs(P_Gen_Total(i)) > 1e-3
                    idList = [idList; i; i];
                    typeList = [typeList; "PG"; "QG"];
                end

                % --- B. Loads ---
                if abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4
                    idList = [idList; i; i];
                    typeList = [typeList; "PD"; "QD"];
                end
            end

            %% 3. Flow Sensors (PF, QF)
            % Logic: Attempt to find a branch for every bus.
            for i = 1:nb
                % A. Try finding an OUTGOING branch (From Bus == i)
                branchIdx = find(obj.mpc.branch(:, F_BUS) == i, 1);

                % B. If no outgoing, try finding an INCOMING branch (To Bus == i)
                if isempty(branchIdx)
                    branchIdx = find(obj.mpc.branch(:, T_BUS) == i, 1);
                end

                % C. Add the Branch ID (Duplicate checks handled at the end)
                if ~isempty(branchIdx)
                    idList = [idList; branchIdx; branchIdx];
                    typeList = [typeList; "PF"; "QF"];
                end
            end

            %% 4. Create, Clean, and Sort Output Table
            sensorTable = table(idList, typeList, 'VariableNames', {'Bus', 'DataType'});

            % CRITICAL FIX: Remove duplicate rows
            % This merges cases where two buses picked the same line.
            sensorTable = unique(sensorTable, 'rows');

            % Sort by ID (mixed Bus and Branch IDs)
            sensorTable = sortrows(sensorTable, 'Bus');
        end

        function sensorTable = createMeasStrategyOld(obj)
            % CREATEMEASSTRATEGY Generates a measurement plan.
            % UPDATED: Excludes PG/QG sensors for Synchronous Condensers (where PG ~ 0)
            % to avoid zero-value issues in GMM training.

            define_constants;
            % mpc = obj.mpc;
            nb = size(obj.mpc.bus, 1);
            nbr = size(obj.mpc.branch, 1);

            idList = [];
            typeList = [];

            % 1. Pre-calculate Total Generation per Bus
            % We need to check the ACTUAL generation value, not just the Bus Type.
            P_Gen_Total = sparse(obj.mpc.gen(:, GEN_BUS), 1, obj.mpc.gen(:, PG), nb, 1);

            %% 1. Voltage Sensors (VM)
            for i = 1:nb
                % Measure VM at all Voltage-Controlled buses (Slack & PV)
                if obj.mpc.bus(i, BUS_TYPE) == 3 || obj.mpc.bus(i, BUS_TYPE) == 2
                    idList = [idList; i];
                    typeList = [typeList; "VM"];
                end
            end

            %% 2. Injection Sensors
            for i = 1:nb
                % --- A. Generators ---
                % Logic: Only measure PG/QG if there is SIGNIFICANT Active Generation.
                % This filters out Synchronous Condensers (Bus 3, 6, 8) where PG is zero.
                if obj.mpc.bus(i, BUS_TYPE) >= 2 && abs(P_Gen_Total(i)) > 1e-3
                    idList = [idList; i; i];
                    typeList = [typeList; "PG"; "QG"];
                end

                % --- B. Loads ---
                % Measure if actual load exists (Filters out empty buses)
                if abs(obj.mpc.bus(i, PD)) > 1e-4 || abs(obj.mpc.bus(i, QD)) > 1e-4
                    idList = [idList; i; i];
                    typeList = [typeList; "PD"; "QD"];
                end

                % Note: We intentionally removed the "Transit Node" (Virtual Zero) logic
                % here to prevent feeding zero values to your GMM.
            end

            %% 3. Flow Sensors (PF, QF)
            % Measure PF/QF on ALL branches to ensure observability.
            % This is critical for Bus 8, which now has NO injection measurements.
            branchIdx = (1:nbr)';

            idList = [idList; branchIdx; branchIdx];
            typeList = [typeList; repmat("PF", nbr, 1); repmat("QF", nbr, 1)];

            %% 4. Create and Sort Output Table
            sensorTable = table(idList, typeList, 'VariableNames', {'Bus', 'DataType'});
            % sensorTable = table(idList, typeList, 'VariableNames', {'ID', 'Type'});

            % Sort by ID Ascending
            sensorTable = sortrows(sensorTable, 'Bus');
            % sensorTable = sortrows(sensorTable, 'ID');
        end
        function sensorTable = createMeasStrategyOldOld(obj)

            % Creates a sensor selection TABLE based on NON-ZERO values in the case file.
            % This table is compatible with the GridEdgeNode.computeMeas() input format.
            %
            % Output Columns:
            %   ID   : Bus Index (Used for ALL types: VM, P/Q Inj, and PF/QF)
            %   Type : String identifier ("VM", "PD", "PF", etc.)

            % 1. Define Constants
            define_constants;
            mpc = obj.mpc;

            % 2. Pre-process Data
            nb = size(mpc.bus, 1);

            % Aggregate Generators per bus
            P_Gen_Total = sparse(mpc.gen(:, GEN_BUS), 1, mpc.gen(:, PG), nb, 1);
            Q_Gen_Total = sparse(mpc.gen(:, GEN_BUS), 1, mpc.gen(:, QG), nb, 1);

            % Find Sending End Buses (Buses that are the "From" end of a line)
            sending_end_buses = unique(mpc.branch(:, F_BUS));

            % Initialize storage vectors
            idList = [];
            typeList = [];

            % 3. Loop through EACH BUS
            for i = 1:nb

                % --- A. Voltage Sensors (VM) ---
                % Logic: Measure V only at Voltage-Controlled Buses (PV or Slack)
                if mpc.bus(i, BUS_TYPE) == 3 || mpc.bus(i, BUS_TYPE) == 2
                    idList = [idList; i];
                    typeList = [typeList; "VM"];
                end

                % --- B. Load Sensors (PD, QD) ---
                % Logic: Place sensor only if actual Load > 1e-4 MW/MVar exists
                if abs(mpc.bus(i, PD)) > 1e-4 || abs(mpc.bus(i, QD)) > 1e-4
                    idList = [idList; i; i];
                    typeList = [typeList; "PD"; "QD"];
                end

                % --- C. Generation Sensors (PG, QG) ---
                % Logic: Place sensor only if Generation > 1e-4 MW/MVar exists
                if abs(P_Gen_Total(i)) > 1e-4 && abs(Q_Gen_Total(i)) > 1e-4
                    idList = [idList; i; i];
                    typeList = [typeList; "PG"; "QG"];
                end

                % --- D. Flow Sensors (PF, QF) ---
                % Logic: If this bus is the start of a line, place a flow sensor.
                % NOTE: This uses the BUS ID as the identifier.
                if ismember(i, sending_end_buses)
                    idList = [idList; i; i];
                    typeList = [typeList; "PF"; "QF"];
                end

            end

            % 4. Create Output Table
            sensorTable = table(idList, typeList, 'VariableNames', {'ID', 'Type'});

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