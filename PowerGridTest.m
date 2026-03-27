classdef PowerGridTest < matlab.unittest.TestCase
    % POWERGRIDTEST Unit tests for the PowerGrid class.
    % Verifies topological observability, measurement strategy generation,
    % and Jacobian rank for different IEEE test cases and strategies.
    
    properties
        Grid14  
        Grid30  
    end
    
    methods (TestMethodSetup)
        function createGridObjects(testCase)
            testCase.Grid14 = PowerGrid('case14');
            testCase.Grid30 = PowerGrid('case30');
        end
    end
    
    methods (Test)
        
        %% =========================================================
        %% 14-BUS TESTS
        %% =========================================================
        function testStrategyGeneration14(testCase)
            testCase.Grid14.createMeasStrategy('case14');
            table = testCase.Grid14.sensorTable;
            testCase.verifyNotEmpty(table, 'Sensor table should not be empty.');
            testCase.verifyTrue(ismember('Bus', table.Properties.VariableNames), ...
                'Sensor table must have a Bus column.');
        end
        
        function testObservability_14bus(testCase)
            grid = testCase.Grid14;
            grid.createMeasStrategy('case14');
            
            H = testCase.buildJacobian(grid);
            H(:,1) = []; % Remove slack bus angle column
            
            nb = size(grid.mpc.bus, 1);
            obsRank = sprank(H);
            reqRank = 2*nb - 1;
            
            fprintf('\n========================================\n');
            fprintf(' IEEE 14 Bus - Original Strategy\n');
            grid.printMeasSummary();
            
            testCase.verifyEqual(obsRank, reqRank, ...
                sprintf('IEEE 14 Jacobian is rank deficient! (Rank: %d, Req: %d)', obsRank, reqRank));
        end

        % =========================================================
        % 30-BUS TESTS: SCENARIO B (RTU HUB)
        % =========================================================
        function testObservability_30bus_Hub(testCase)
            grid = testCase.Grid30;
            grid.createMeasStrategy('case30');
            
            H = testCase.buildJacobian(grid);
            H(:,1) = []; 
            
            nb = size(grid.mpc.bus, 1);
            obsRank = sprank(H);
            reqRank = 2*nb - 1;
            
            fprintf('\n========================================\n');
            fprintf(' IEEE 30 Bus - RTU HUB Strategy\n');
            grid.printMeasSummary();
            
            testCase.verifyEqual(obsRank, reqRank, ...
                sprintf('IEEE 30 Hub Jacobian is rank deficient! (Rank: %d, Req: %d)', obsRank, reqRank));
        end
        
        % =========================================================
        % DATA EXTRACTION COMPATIBILITY TEST
        % =========================================================
        function testDataExtractionMapping30(testCase)
            % Verifies that collectBusData successfully builds columns for
            % EVERY measurement assigned by the strategy, ensuring no crashes.
            
            strategies = {'thesis', 'hub'};
            
            for s = 1:length(strategies)
                strat = strategies{s};
                grid = PowerGrid('case30');
                grid.createMeasStrategy('case30');
                
                % 1. Run power flow to get the true results
                mpopt = mpoption('out.all', 0, 'verbose', 0);
                PowerFlowResults = runpf(grid.mpc, mpopt);
                
                % 2. Collect the true data using F_BUS strictly
                busData = grid.collectBusData(PowerFlowResults);
                
                % 3. Extract the measurement strategy table
                tbl = grid.sensorTable;
                flowSensors = tbl(tbl.DataType == "PF" | tbl.DataType == "QF", :);
                
                % 4. Verify all flow columns exist and check for 0-variance ML risk
                for k = 1:height(flowSensors)
                    busID = flowSensors.Bus(k);
                    mType = string(flowSensors.DataType(k));
                    brID = flowSensors.Branch(k);
                    
                    colName = sprintf('%s-%d', mType, brID);
                    
                    try
                        % If collectBusData failed to create the column, this will crash
                        actualValue = busData(busID).trueVal{1, colName};
                    catch
                        testCase.verifyFail(sprintf('[%s Strategy] CRASH: Missing column %s at Bus %d. F_BUS mapping failed.', upper(strat), colName, busID));
                        continue;
                    end
                    
                    % Soft Warning: If the flow is exactly 0, alert the user about the ML noise floor
                    if actualValue == 0
                        fprintf('[ML WARNING] %s measured at Bus %d in %s strategy is exactly 0. Ensure noise floor is applied in addThermalNoise().\n', colName, busID, upper(strat));
                    end
                end
            end
        end
    end
    
    % =====================================================================
    % Helper Methods
    % =====================================================================
    methods (Access = private)
        function H = buildJacobian(~, grid)
            % BUILDJACOBIAN Constructs a structural flat-start Jacobian 
            % to test for topological observability.
            
            define_constants;
            mpc = grid.mpc;
            nb = size(mpc.bus, 1);
            
            % Generate Admittance Matrices using MATPOWER
            [Ybus, Yf, Yt] = makeYbus(mpc);
            
            % Extract Conductance (G) and Susceptance (B)
            G = real(Ybus); B = imag(Ybus);
            Gf = real(Yf);  Bf = imag(Yf);
            Gt = real(Yt);  Bt = imag(Yt);
            
            tbl = grid.sensorTable;
            nm = height(tbl);
            
            % Initialize H matrix: [dZ/dTheta, dZ/dV]
            H = zeros(nm, 2*nb); 
            
            for k = 1:nm
                type = string(tbl.DataType(k));
                bus = tbl.Bus(k);
                branch = tbl.Branch(k);
                
                if type == "VM"
                    % dV_i / dV_i = 1
                    H(k, nb + bus) = 1;
                    
                elseif type == "PG" || type == "PD"
                    % Active Power Injection
                    H(k, 1:nb) = -B(bus, :);      % dP / dTheta
                    H(k, nb+1:end) = G(bus, :);   % dP / dV
                    
                elseif type == "QG" || type == "QD"
                    % Reactive Power Injection
                    H(k, 1:nb) = G(bus, :);       % dQ / dTheta
                    H(k, nb+1:end) = -B(bus, :);  % dQ / dV
                    
                elseif type == "PF"
                    % Active Power Flow
                    if bus == mpc.branch(branch, F_BUS)
                        H(k, 1:nb) = -Bf(branch, :);
                        H(k, nb+1:end) = Gf(branch, :);
                    else
                        H(k, 1:nb) = -Bt(branch, :);
                        H(k, nb+1:end) = Gt(branch, :);
                    end
                    
                elseif type == "QF"
                    % Reactive Power Flow
                    if bus == mpc.branch(branch, F_BUS)
                        H(k, 1:nb) = Gf(branch, :);
                        H(k, nb+1:end) = -Bf(branch, :);
                    else
                        H(k, 1:nb) = Gt(branch, :);
                        H(k, nb+1:end) = -Bt(branch, :);
                    end
                end
            end
        end
    end
end