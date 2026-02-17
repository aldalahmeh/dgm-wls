classdef TestGridEdgeNode < matlab.unittest.TestCase
    %TESTGRIDEDGENODE Unit test for State Estimation measurement function h(x)
    %   Verifies that computeMeas returns correct Net Injections and Line Flows
    %   matching MATPOWER's Ground Truth.

    properties
        Grid        % PowerGrid Object
        EdgeNode    % GridEdgeNode Object
        x_true      % Table containing True State (VM, VA)
        nBus        % Number of buses
        baseMVA     % System Base MVA (for unit conversion)
    end

    methods(TestMethodSetup)
        function setupEnvironment(testCase)
            % 1. Load the Case (IEEE 14)
            casename = 'case14';
            testCase.Grid = PowerGrid(casename);
            testCase.nBus = numel(testCase.Grid.Bus);

            % 2. Initialize Edge Node
            % (nSim=1 is sufficient for physics validation)
            testCase.EdgeNode = GridEdgeNode(testCase.Grid, 1);

            % 3. Extract True State
            % PowerGrid automatically runs runpf on initialization
            % and stores results in TrueGridState (VM is p.u., VA is Degrees)
            testCase.x_true = testCase.Grid.TrueGridState;
            testCase.baseMVA = testCase.Grid.mpc.baseMVA;

            % Verify x_true has the correct columns for computeMeas
            % computeMeas expects table with 'VM' and 'VA'
            assert(ismember('VM', testCase.x_true.Properties.VariableNames), 'x_true missing VM');
            assert(ismember('VA', testCase.x_true.Properties.VariableNames), 'x_true missing VA');
        end
    end

    methods(Test)

        function testNetInjections(testCase)
            % TESTNETINJECTIONS Verify P and Q Net Injections at all buses
            % h(x) for 'PG'/'PD' should return Net P Injection
            % h(x) for 'QG'/'QD' should return Net Q Injection
            % h(x) for 'PF' should return flow P Injection
            % h(x) for 'QF' should return flow Q Injection

            % Construct Input z (Measurement Request)
            zRequest = testCase.Grid.createMeasStrategy();

            % Run h(x)
            zHat = testCase.EdgeNode.computeMeas(testCase.x_true, zRequest);
            zHatTable = testCase.meas2busTable(zHat);

            % Compute Ground Truth (Net Injection)
            % Use zRequest to extract required ground truth
            expected_P = testCase.extractExptdVals(zRequest,'PD');
            expected_Q = testCase.extractExptdVals(zRequest,'QD');

            % Compute Ground Truth (Fllow Injection)
            % Use zRequest to extract required ground truth
            expected_PF = testCase.extractExptdVals(zRequest,'PF');
            expected_QF = testCase.extractExptdVals(zRequest,'QF');


            % Verify
            actual_P = zHatTable.PD;
            actual_Q = zHatTable.QD;
            actual_PF = zHatTable.PF;
            actual_QF = zHatTable.QF;

            % Check PD
            testCase.verifyEqual(actual_P, expected_P, 'AbsTol', 1e-4, ...
                'Active Power Injection (P) mismatch');

            % Check QD
            testCase.verifyEqual(actual_Q, expected_Q, 'AbsTol', 1e-4, ...
                'Reactive Power Injection (Q) mismatch');

            % Check PF
            testCase.verifyEqual(actual_PF, expected_PF, 'AbsTol', 1e-4, ...
                'Active Power Flow (PF) mismatch');

            % Check QF
            testCase.verifyEqual(actual_QF, expected_QF, 'AbsTol', 1e-4, ...
                'Reactive Power Flow (QF) mismatch');
        end

        function testJacbianMatrix(testCase)


            nBus = numel(testCase.Grid.Bus);

            % Construct Input z (Measurement Request)
            zRequest = testCase.Grid.createMeasStrategy();

            % Flat start
            x_flat = testCase.x_true;
            x_flat.VA(:) = 0;
            x_flat.VM(:) = 1.0;

            % Compute the Jacobian matrix using the EdgeNode's method
            H = testCase.EdgeNode.computeJacobian(x_flat, zRequest);


            % 3. Check Dimensions
            testCase.verifySize(full(H), [height(zRequest),2*nBus], 'H matrix has wrong size');

            % 4. Verify Values
            % The first nBus columns (Angles) must be EXACTLY 0
            testCase.verifyEqual(full(H(1, 1:nBus)), zeros(1, nBus), ...
                'VM derivative w.r.t Angle must be 0');

            % The column for V_1 (index nBus + 1) must be EXACTLY 1
            testCase.verifyEqual(full(H(1, nBus+1)), 1, ...
                'VM derivative w.r.t Self Voltage must be 1');

            % 5. Check the rank
            % Rank should be #Buses - 2
            testCase.verifyGreaterThanOrEqual(rank(full(H)), 2*nBus - 1, ...
                'Rank deficient! rank(H) >= 2*#Buses - 1');



        end

        function testWlsConvergenceFlatStart(testCase)
            % testWlsConvergenceFlatStart Verify WLS recovers True State from Flat Start
            % using Perfect Measurements (No Noise).

            fprintf('\n--- Running WLS Flat Start Test ---\n');

            % 1. Get Ground Truth (xTrue)
            % This is the exact solution from Power Flow (runpf)
            xTrue = testCase.x_true;

            % 2. Generate Perfect Measurements (No Noise)
            % We use the strategy, but calculate z based on xTrue
            measStrategy = testCase.Grid.createMeasStrategy();
            zPerfect = testCase.EdgeNode.computeMeas(xTrue, measStrategy);

            % Artificial Zero Variance for the test logic (High Weight)
            % We simulate "Perfect Sensors" by giving them tiny variance in the table
            % (Note: In the WLS code, W = 1/Variance)
            zPerfect.Variances = 1e-6 * ones(height(zPerfect), 1);

            % No corrupted data
            zPerfect.Corrupt = zeros(height(zPerfect),1);

            % 3. Set Initial Estimate to Flat Start
            % (Reset the object's internal estimate to V=1, Angle=0)
            % nBus = height(xTrue);
            % testCase.EdgeNode.initEst = table(ones(nBus,1), zeros(nBus,1), ...
            %     'VariableNames', {'VM', 'VA'});

            % 4. Run WLS with Perfect Data
            % We pass zPerfect to override the internal noisy generation
            xEst = testCase.EdgeNode.computeStateEst(zPerfect);

            % 5. Validation
            % A) Check Voltage Magnitudes (Target: < 0.001 p.u. error)
            errorVm = max(abs(xEst.VM - xTrue.VM));
            fprintf('Max Voltage Error: %.6f p.u.\n', errorVm);
            testCase.verifyLessThan(errorVm, 1e-3, 'Voltage Estimation failed to converge');

            % B) Check Voltage Angles (Target: < 0.1 deg error)
            % Note: Angles are relative. If the slack bus reference shifted,
            % we compare the *differences* relative to slack, or just absolute if slack is fixed.
            % Since we fixed Slack=0 in WLS, absolute comparison is valid.
            errorVa = max(abs(xEst.VA - xTrue.VA));
            fprintf('Max Angle Error:   %.6f deg\n', errorVa);
            testCase.verifyLessThan(errorVa, 0.1, 'Angle Estimation failed to converge');

            % 6. Print Comparison Table
            fprintf('\n%5s | %15s | %15s || %15s | %15s\n', ...
                'Bus', 'VM True', 'VM Est', 'VA True', 'VA Est');
            fprintf('%s\n', repmat('-', 1, 75));

            for i = 1:nBus
                fprintf('%5d | %15.4f | %15.4f || %15.4f | %15.4f\n', ...
                    i, ...
                    xTrue.VM(i), xEst.VM(i), ...
                    xTrue.VA(i), xEst.VA(i));
            end
            fprintf('%s\n', repmat('-', 1, 75));
        end

        function testInitEst_Flat(testCase)
            % TESTINITEST_FLAT Verify 'flat' string returns V=1, Ang=0

            % 1. Call method
            % (Measurements are ignored for flat start, pass empty or random)
            initEst = testCase.EdgeNode.computeInitEst('flat');

            % 2. Verify
            testCase.verifyEqual(initEst.VM, ones(testCase.nBus, 1), ...
                'Flat start VM must be all 1.0');
            testCase.verifyEqual(initEst.VA, zeros(testCase.nBus, 1), ...
                'Flat start VA must be all 0.0');
        end
        % 
        % function testInitEstWarmZeroInjection(testCase)
        %     % TESTINITEST_WARM_ZEROINJECTION
        %     % Case: Zero P and Q load everywhere.
        %     % Expectation: Result should effectively be Flat Start.
        % 
        %     % 1. Create a "Zero" Measurement Set
        %    zMeas = testCase.createMockMeas();
        % 
        %     % 2. Run Warm Start
        %     initEst = testCase.EdgeNode.computeInitEst('warm', zMeas);
        % 
        %     % 3. Verify
        %     % Angles should be 0 (allow tiny tolerance for matrix inversion noise)
        %     testCase.verifyEqual(initEst.VA, zeros(testCase.nBus, 1), 'AbsTol', 1e-6, ...
        %         'Zero P injection should result in 0 angles');
        % 
        %     % Voltages should be 1 (Note: PV buses might be set to V_setpoint in your logic,
        %     % but if Q=0 and V_set=1, it stays 1. Check PQ buses mainly.)
        %     testCase.verifyEqual(initEst.VM, ones(testCase.nBus, 1), 'AbsTol', 1e-6, ...
        %         'Zero Q injection should result in flat voltages');
        % end
        % 
        % function testInitEstWarmPhysicsP(testCase)
        %     % TESTINITEST_WARM_PHYSICSP
        %     % Case: Inject Active Load (PD) at Bus 2. Zero Q.
        %     % Expectation: Bus 2 Angle should go NEGATIVE. Voltages constant.
        % 
        %     % 1. Create Measurement Set
        %     zMeas = testCase.createMockMeas();
        % 
        %     % 2. Add Active Load at Bus 2 (e.g., 1.0 p.u. MW)
        %     % Find row for Bus 2, Type PD
        %     idxPD = find(zMeas.Bus == 2 & strcmp(zMeas.DataType, 'PD'));
        %     if ~isempty(idxPD)
        %         zMeas.Measurements(idxPD) = 1.0; % Positive Load = Negative Injection
        %     end
        % 
        %     % 3. Run Warm Start
        %     initEst = testCase.EdgeNode.computeInitEst('warm', zMeas);
        % 
        %     % 4. Verify Decoupling
        %     % Angle at Bus 2 should be NEGATIVE (Power flows 1 -> 2)
        %     angleBus2 = initEst.VA(2);
        %     testCase.verifyLessThan(angleBus2, -0.1, ...
        %         'Active Load should cause negative voltage angle');
        % 
        %     % Voltages should generally stay 1.0 (Pure DC flow assumption)
        %     % (Allowing slack for PV bus setpoints if they differ)
        %     testCase.verifyEqual(initEst.VM(2), 1.0, 'AbsTol', 1e-4, ...
        %         'Pure Active Power load should not change Voltage Magnitude in decoupled init');
        % end

    end

    % Helper functions
    methods
        function zMeas = createMockMeas(obj)

            zMeas = obj.Grid.createMeasStrategy();
            zMeas = renamevars(zMeas,'ID','Bus');
            zMeas = renamevars(zMeas,'Type','DataType');
            zMeas.Measurements = zeros(height(zMeas), 1); % FORCE ALL TO 0
            zMeas.Corrupt = zeros(height(zMeas), 1);

        end
        function zHatTable = meas2busTable(obj, zHat)


            % Preallocate table
            nBus = numel(obj.Grid.Bus);
            baseVal = obj.Grid.baseVal;

            Bus = (1:nBus).';
            VM = nan(nBus,1);
            VA = nan(nBus,1);
            PD = nan(nBus,1);
            PQ = nan(nBus,1);
            PG = nan(nBus,1);
            QG = nan(nBus,1);
            PF = nan(nBus,1);
            QF = nan(nBus,1);

            zHatTable = table(Bus, VA, VM, PD, PQ, PG, QG, PF, QF, ...
                'VariableNames', {'Bus', 'VA', 'VM', 'PD','QD','PG','QG','PF','QF'});

            % Polpulate table
            nMeas = numel(zHat(:,1));

            for iMeas = 1:nMeas
                busNum = table2array(zHat(iMeas,1));
                measType = table2array(zHat(iMeas,2));
                zHatTable{busNum, measType} = table2array(zHat(iMeas, 3))*baseVal;
            end

        end

        function exptcdVal = extractExptdVals(obj, Request, measType)

            % Get all relevant data
            switch measType
                case 'PD'
                    exptcdData = obj.Grid.TrueNetPower.PD;
                case 'QD'
                    exptcdData = obj.Grid.TrueNetPower.QD;
                case 'PF'
                    exptcdData = obj.Grid.TruePowerFlow.PF;
                case 'QF'
                    exptcdData = obj.Grid.TruePowerFlow.QF;

            end

            % Extract data as per the Request input
            % Reduce the Request to nBus logical array to slice the
            % required data

            nBus = numel(obj.Grid.Bus);
            reqDataInd = zeros(1,nBus);
            exptcdVal = nan*ones(1,nBus);

            for iBus = 1:nBus
                % Check if the measType is there in iBus
                temp = Request{(Request{:,1} == iBus),2};
                reqDataInd(iBus) = any(temp == measType);

            end
            reqDataInd = logical(reqDataInd);
            exptcdVal(reqDataInd) = exptcdData(reqDataInd);
            exptcdVal = exptcdVal';

        end

    end
end