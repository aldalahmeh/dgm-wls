%% Notes
% 1- Collect sensor measurements. (DONE)
% 2- Get sensor nodes variances. (DONE)
% 3- Compute initial estimate. (DONE)
% 4- Compute WLS estimate. (DONE)

%%
classdef GridEdgeNode < handle
    % GRIDEDGENODE Summary of this class goes here
    %   Detailed explanation goes here

    properties
        nSim
        Grid
        sensorData
        sensorMeas
        wlsParam
        DEBUG
    end

    methods
        function obj = GridEdgeNode(Grid, nSim, maxIter, tol)
            %GridEdgeNode Construct an instance of this class
            %   Detailed explanation goes here

            % Debug flag
            obj.DEBUG = false;

            obj.nSim = nSim;
            obj.Grid = Grid;

            % WLS parameters
            if nargin == 4
                obj.wlsParam.maxIter = maxIter;
                obj.wlsParam.tol = tol;
            elseif  nargin == 3
                obj.wlsParam.maxIter = maxIter;
                obj.wlsParam.tol = 1e-4;
            elseif nargin < 3
                obj.wlsParam.maxIter = 20;
                obj.wlsParam.tol = 1e-4;
            end

            obj.wlsParam.alpha = 1;

            % Prepare the Bp matrix
            obj.prepareDecoupledMat();

            nBus = numel(Grid.Bus);

            % Preallocate sensorData struct
            obj.sensorData = repmat(struct('Bus_ID', [], 'Chains', []), 1, nBus );

            % Define the column names as a cell array or string array
            varNames = {'Bus #', 'Data', 'Data-Type'};

            % Create the empty table
            obj.sensorMeas = table([],[],[],'VariableNames', varNames);

        end

        function receive(obj, rxData)

            % Get sensor ID
            sensorID = rxData.Bus_ID;

            obj.sensorData(sensorID) = rxData;

        end

        function initEst = computeInitEst(obj, initType, measurements )
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            nBus = numel(obj.Grid.Bus);

            if strcmpi(initType,'flat')

                vEstArray = ones(nBus,1);
                vAngEstArray = zeros(nBus, 1);

            elseif strcmpi(initType,'warm')


                % Get power injection measurements
                [Pinj, Qinj] = obj.collectPQinj(measurements);

                % Get admittance derivative matrices
                Bp  = obj.Grid.Bp;
                Bpp = obj.Grid.Bpp;

                % Compute linreazied initial theta
                vAngEst = rad2deg(-Bp \ Pinj);

                % Append a zero for the reference (slack) bus value
                nBus = numel(obj.Grid.Bus);
                slackBusInd = obj.getBusTypeInd('Ref');
                nonSlackIdx = setdiff(1:nBus, slackBusInd);
                vAngEstArray = zeros(nBus, 1);
                vAngEstArray(nonSlackIdx) = vAngEst;

                % Compute the linearized initial voltage
                vEst  = 1.0 + Bpp \ Qinj;
                vEstArray = ones(nBus,1);
                pqBusInd = obj.getBusTypeInd('PQ');
                vEstArray(pqBusInd) = vEst;


            else
                error('Undefined state initialization type!')

            end

            % Store in a table
            initEst = table((1:nBus)', vEstArray, vAngEstArray,...
                'VariableNames',{'Bus_ID', 'VM', 'VA'});
        end


        function [stateEst, convergenceData] = computeStateEst(obj, Options)
            % COMPUTESTATEEST Solves WLS/IRWLS State Estimation using Gauss-Newton
            
            % Collect parameters
            busMeasMapping = obj.Grid.sensorTable;
            maxIter = obj.wlsParam.maxIter;
            tol     = obj.wlsParam.tol;
            alpha   = obj.wlsParam.alpha;
            nSim    = obj.nSim;

            % Initialize output
            stateEst     = cell(nSim,1);
            convergenceData  = nan(nSim,1);

            
            % Setup Parallel / Debug Mode and DataQueue
            q = parallel.pool.DataQueue;
            afterEach(q, @(data) obj.updateMCProgress(data, nSim));
            
            if isfield(Options, 'DebugMode') && Options.DebugMode == true
                numWorkers = 0; % Forces sequential execution
                fprintf('--- RUNNING IN DEBUG MODE (Sequential) ---\n');
            else
                numWorkers = Inf; % Uses the full parallel pool
                fprintf('Starting Monte Carlo simulation (%d runs)...\n', nSim);
            end
           
            % Parallel loop
            parfor (iSim = 1:nSim, numWorkers)
            % for iSim = 1:nSim
                % Silence parfor warnings   
                wlsConverged = false;       
                xCurrent = []; 
                
                % Prepare measurements
                zTable = obj.prepareSensorMeas(busMeasMapping, iSim);
                zTable = obj.filterOverride(zTable, Options.FilteringOverride);

                % Estimate corrupted data if enabled
                if ~strcmp(Options.dataImpute,'None') && any(zTable.Corrupt)
                    zTable = obj.estCorruptData(zTable, Options);
                end
                
                % Compute initial estimate
                if strcmp(Options.WlsInit, 'warm')
                    xCurrent = obj.computeInitEst('Warm', zTable);
                elseif strcmp(Options.WlsInit, 'flat')
                    xCurrent = obj.computeInitEst('Flat', zTable);
                else 
                    error('Unknown initialization method!')
                end
                
                % Get comms noise variances and compute BASE weight matrix                
                commVarTable = obj.prepareCommVar(iSim);
                W_base = obj.computeWeightMatrix(zTable, commVarTable, Options);
                
                % Number of buses
                nBus = height(xCurrent);
                
                % 2. Iteration Loop
                for iter = 1:maxIter
                    % Compute zHat = h(x)
                    zHatTable = obj.computeMeas(xCurrent, busMeasMapping);
                    zHatTable.Corrupt = zTable.Corrupt;
                    zHatTable.Override = zTable.Override;
                    
                    % Compute residuals
                    r = obj.computeResiduals(zTable, zHatTable);
                  
                    % --- NEW: DYNAMIC WEIGHT UPDATE FOR IRWLS ---
                    % If Options.RobustEstimator = 'Huber', this updates the weights.
                    % Otherwise, it just returns W_base for ordinary WLS.
                    W_iter = obj.applyRobustWeights(W_base, r.Residuals, Options);
                    
                    % Construct the Jacobian matrix
                    H = obj.computeJacobian(xCurrent, zHatTable);
                    
                    % Remove slack bus angle column
                    Hred = H;
                    Hred(:,1) = [];
                    
                    % Solve normal equations using the current iteration's weights
                    G = Hred' * W_iter * Hred;
                    rhs = Hred' * W_iter * r.Residuals;
                    
                    % Solve for reduced update vector
                    lastwarn('');
                    dx_red = G \ rhs;
                  
                    [~, msgid] = lastwarn;
                    if strcmp(msgid, 'MATLAB:nearlySingularMatrix')
                        error('WLS Failed: Gain matrix is singular. Check Observability (Rank < 2*nBus-1).');
                    end
                    
                    % Reconstruct Full Update Vector
                    dx = alpha * [0; dx_red];
                    
                    % Update State 
                    dVa_deg = rad2deg(dx(1:nBus));      
                    dVm     = dx(nBus+1:end);  

                    % --- NEW: STEP LIMITER (DAMPING) ---
                    % Prevent massive non-linear overshoots caused by impulsive noise
                    % max_dVa = 15.0; % Maximum angle change per iteration (Degrees)
                    % max_dVm = 0.15; % Maximum voltage magnitude change per iteration (p.u.)
                    % 
                    % dVa_deg = max(min(dVa_deg, max_dVa), -max_dVa);
                    % dVm     = max(min(dVm, max_dVm), -max_dVm);
                    
                    xCurrent.VA = xCurrent.VA + dVa_deg;
                    xCurrent.VM = xCurrent.VM + dVm;
                    
                    % Check Convergence
                    max_update = max(abs(dx));

                    if max_update < tol
                        wlsConverged = true;
                        % Store # convergence iterations
                        convergenceData(iSim) = iter;

                        % fprintf('WLS Converged in %d iterations.\n', iter);
                        stateEst{iSim} = xCurrent;
                        
                        % Only print in debug mode to keep parfor fast
                        if numWorkers == 0
                            fprintf('WLS Converged in %d iterations (Sim %d).\n', iter, iSim);
                        end
                        break
                    end
                end
                
                if ~wlsConverged
                    warnMsg = sprintf('WLS did not converge within %d iterations for run %d.', maxIter, iSim);
                    send(q, warnMsg);                    
                    stateEst{iSim} = xCurrent; % Return best effort
                else
                    send(q, iSim); % Send progress update
                end
            end
        end

        function stateEst = computeStateEstOld(obj, Options)
            % COMPUTESTATEEST Solves WLS State Estimation using Gauss-Newton method

            % Measurements-bus mapping
            busMeasMapping = obj.Grid.sensorTable;
            maxIter = obj.wlsParam.maxIter;
            tol     = obj.wlsParam.tol;
            alpha   = obj.wlsParam.alpha;

            % Get number of data points
            nSim = obj.nSim;

            % Initialize state estimation
            stateEst = cell(nSim,1);

            % 1. Setup the Parallel Data Queue for Logging
            % fprintf('Starting Monte Carlo simulation (%d runs)...\n', nSim);
            q = parallel.pool.DataQueue;

            % Define what happens when a worker sends a message back.
            % We will use a custom local function (defined at the bottom) to print progress.
            afterEach(q, @(data) obj.updateMCProgress(data, nSim));

            % Determine parallel execution mode
            if isfield(Options, 'DebugMode') && Options.DebugMode == true
                numWorkers = 0; % Forces sequential execution (allows breakpoints)
                fprintf('--- RUNNING IN DEBUG MODE (Sequential) ---\n');
            else
                numWorkers = Inf; % Uses the full parallel pool
                fprintf('Starting Monte Carlo simulation (%d runs)...\n', nSim);                
            end

          
            % Parrel loop
            % parfor (iSim = 1:nSim, numWorkers)
            for iSim = 1:nSim

                % Convergenc flag   
                wlsConverged = false;       

                xCurrent = []; % Initialize to satisfy the static analyzer

                % Print iteration number
                fprintf('MC iter #: %d\n',iSim);

                % Prepare measurements
                zTable = obj.prepareSensorMeas(busMeasMapping, iSim);

                % Check if an override is needed
                zTable = obj.filterOverride(zTable, Options.FilteringOverride);

                % Compute initial estimate
                if strcmp(Options.WlsInit, 'warm')
                    xCurrent = obj.computeInitEst('Warm', zTable);
                elseif strcmp(Options.WlsInit, 'flat')
                    xCurrent = obj.computeInitEst('Flat', zTable);
                else 
                    error('Unkown initilization method!')
                end

                % Get comms noise variances
                commVarTable = obj.prepareCommVar(iSim);

                % Compute weight matrix                
                W = obj.computeWeightMatrix(zTable, commVarTable, Options);
                % W = inv(diag(zTable.Variances) + diag(commVarTable.Comm_Variance));

                
                % WLS iteration
                % Number of buses
                nBus = height(xCurrent);

                % 2. Iteration Loop
                fprintf('Starting WLS Estimation...\n');
               

                for iter = 1:maxIter
                    % Compute zHat = h(x)
                    zHatTable = obj.computeMeas(xCurrent, busMeasMapping);
                    zHatTable.Corrupt = zTable.Corrupt;
                    zHatTable.Override = zTable.Override;

                    % Compute residuals
                    r = obj.computeResiduals(zTable, zHatTable);
                  
                    % Construct the Jacobian matrix
                    H = obj.computeJacobian(xCurrent, zHatTable);

                    % Remove slack bus angle column
                    Hred = H;
                    Hred(:,1) = [];

                    % Solve normal equations
                    G = Hred' * W * Hred;

                    % RHS = H_red' * W * r
                    rhs = Hred' * W * r.Residuals;

                    % Solve for reduced update vector
                    % Use decomposition or direct solve. Warning suppression for singular cases.
                    lastwarn('');
                    dx_red = G \ rhs;
                  
                    [~, msgid] = lastwarn;
                    if strcmp(msgid, 'MATLAB:nearlySingularMatrix')
                        error('WLS Failed: Gain matrix is singular. Check Observability (Rank < 2*nBus-1).');
                    end

                    % F. Reconstruct Full Update Vector
                    % Insert 0 for the Slack Angle update
                    dx = alpha * [0; dx_red];

                    % G. Update State (With Unit Conversion)
                    % H is usually computed w.r.t RADIANS (MATPOWER standard).
                    % x.VA is usually in DEGREES.
                    % We must convert dVa from Radians to Degrees.

                    dVa_rad = dx(1:nBus);      % First nBus elements are Angles
                    dVm     = dx(nBus+1:end);  % Remaining elements are Voltages

                    dVa_deg = rad2deg(dVa_rad);

                    % Apply Update
                    xCurrent.VA = xCurrent.VA + dVa_deg;
                    xCurrent.VM = xCurrent.VM + dVm;

                    % H. Check Convergence
                    max_update = max(abs(dx));
                    %fprintf('Iter %d: Max Update = %e, Cost J = %f\n', iter, max_update, J);

                    if max_update < tol
                        wlsConverged = true;
                        fprintf('WLS Converged in %d iterations.\n', iter);
                        stateEst{iSim} = xCurrent;

                        % Exit the WLS loop
                        break
                    end
                end

                if ~wlsConverged
                    warnMsg = sprintf('WLS did not converge within %d iterations for run %d.', maxIter, iSim);
                    send(q, warnMsg);    
                    stateEst{iSim} = xCurrent; % Return best effort
                end
            end
        end

        % Helper functions
        function estData = estCorruptData(obj, data, Options)
            % ESTCORRUPTDATA Replaces corrupted data with a synthesized estimate.
            % Selects between 'Historical' (Phantom Anchor) and 'MAP' (Hybrid Fusion)
            % based on the Options.dataImpute parameter.
            
            estData = data;
            basVal = obj.Grid.baseVal; 
            
            % 1. Get the actual ROW INDICES of the corrupted measurements
            corruptRows = find(data.Corrupt);
            
            for k = 1:length(corruptRows)
                rowIdx = corruptRows(k);
                
                busID = data.Bus(rowIdx);
                dataType = string(data.DataType(rowIdx)); 
                
                % Skip PG/QG if they were already merged into PD/QD
                if dataType == "PG" || dataType == "QG"
                    continue;
                end
                
                % Get the sensor chains for this specific bus
                sensorChains = obj.sensorData(busID).Chains;
                if isempty(sensorChains)
                    continue;
                end
                
                % Map the DataType to the correct Chain Name and Dimension (1=P/V, 2=Q)
                if dataType == "VM"
                    chainName = "VM";
                    dimIdx = 1;
                elseif dataType == "PD"
                    chainName = "PD-QD";
                    dimIdx = 1;
                elseif dataType == "QD"
                    chainName = "PD-QD";
                    dimIdx = 2;
                elseif startsWith(dataType, "PF")
                    chainName = "PF-QF-" + extractAfter(dataType, "-");
                    dimIdx = 1;
                elseif startsWith(dataType, "QF")
                    chainName = "PF-QF-" + extractAfter(dataType, "-");
                    dimIdx = 2;
                else
                    continue; % Unhandled type
                end
                
                % Find the corresponding chain object
                chainIdx = find([sensorChains.dataType] == chainName);
                if isempty(chainIdx)
                    continue; 
                end
                chainObj = sensorChains(chainIdx);
                
                % Extract GMM parameters
                mu     = chainObj.Means;                
                prob   = chainObj.prob;
                sigma2 = chainObj.Variances;
                
                % 2. Identify Clean and Noise Components
                [~, cleanIdx] = max(prob); % The healthy historical cluster
                [~, noiseIdx] = min(prob); % The cluster representing the corruption
                
                % Extract parameters for the specific dimension
                if chainName == "VM"
                    mu_clean  = mu(cleanIdx, 1);
                    var_clean = sigma2(1, 1, cleanIdx);
                    var_noise = sigma2(1, 1, noiseIdx);
                else
                    mu_clean  = mu(cleanIdx, dimIdx);
                    var_clean = sigma2(dimIdx, dimIdx, cleanIdx);
                    var_noise = sigma2(dimIdx, dimIdx, noiseIdx);
                end
                
                % ---------------------------------------------------------
                % ROUTING: Select Imputation Methodology
                % ---------------------------------------------------------
                if strcmpi(Options.dataImpute, 'Historical')
                    % --- ITEM 4: The "Phantom Anchor" Method ---
                    if chainName == "VM"
                        % Voltages are stable; use the tight clean variance
                        estData.Measurements(rowIdx) = mu_clean;
                        estData.Variances(rowIdx)    = var_clean;
                    else
                        % Power is dynamic; use clean mean but massive noise variance
                        zPseudo = mu_clean / basVal;
                        if dataType == "PD" || dataType == "QD"
                            zPseudo = -abs(zPseudo); 
                        end
                        estData.Measurements(rowIdx) = zPseudo;
                        
                        % Mathematically mimics Huber down-weighting
                        varPseudo = var_noise / (basVal^2);
                        estData.Variances(rowIdx) = varPseudo; 
                    end
                    
                elseif strcmpi(Options.dataImpute, 'MAP')
                    % --- ITEM 3: Hybrid MAP-Style Fusion ---
                    
                    % Un-scale the real-time measurement back to physical units
                    z_meas_raw = data.Measurements(rowIdx);
                    if chainName ~= "VM"
                        z_meas_raw = z_meas_raw * basVal;
                        if dataType == "PD" || dataType == "QD"
                            z_meas_raw = abs(z_meas_raw); 
                        end
                    end
                    
                    % The Weighted Average Value (Pulling wild measurement to safe mean)
                    weight_clean = var_noise / (var_clean + var_noise);
                    weight_noise = var_clean / (var_clean + var_noise);
                    zPseudo = (weight_clean * mu_clean) + (weight_noise * z_meas_raw);
                    
                    % The "Low Confidence" Variance (Elastic band)
                    varPseudo = var_noise; 
                    
                    % Re-apply WLS Scaling and Sign Conventions
                    if chainName ~= "VM"
                        zPseudo = zPseudo / basVal;
                        varPseudo = varPseudo / (basVal^2);
                        if dataType == "PD" || dataType == "QD"
                            zPseudo = -abs(zPseudo); 
                        end
                    end
                    
                    estData.Measurements(rowIdx) = zPseudo;
                    estData.Variances(rowIdx)    = max(varPseudo, 1e-6);
                    
                else
                    error('Unknown imputation method. Set Options.dataImpute to ''Historical'' or ''MAP''.');
                end
                
                % 7. CRITICAL: Un-flag the row so the WLS solver keeps it in the Jacobian!
                % (Commented out as the weight matrix function presently handles this execution).
                % estData.Corrupt(rowIdx) = false;
            end
        end

        function estData = estCorruptData2(obj, data)
            % ESTCORRUPTDATA Replaces corrupted data with a hybrid MAP-style estimate
            % using GMM historical means and real-time noisy measurements.
            
            estData = data;
            basVal = obj.Grid.baseVal; 
            
            % 1. Get the actual ROW INDICES of the corrupted measurements
            corruptRows = find(data.Corrupt);
            
            for k = 1:length(corruptRows)
                rowIdx = corruptRows(k);
                
                busID = data.Bus(rowIdx);
                dataType = string(data.DataType(rowIdx)); 
                
                % Skip PG/QG if they were already merged into PD/QD
                if dataType == "PG" || dataType == "QG"
                    continue;
                end
                
                % Get the sensor chains for this specific bus
                sensorChains = obj.sensorData(busID).Chains;
                if isempty(sensorChains)
                    continue;
                end
                
                % Map the DataType to the correct Chain Name and Dimension (1=P/V, 2=Q)
                if dataType == "VM"
                    chainName = "VM";
                    dimIdx = 1;
                elseif dataType == "PD"
                    chainName = "PD-QD";
                    dimIdx = 1;
                elseif dataType == "QD"
                    chainName = "PD-QD";
                    dimIdx = 2;
                elseif startsWith(dataType, "PF")
                    chainName = "PF-QF-" + extractAfter(dataType, "-");
                    dimIdx = 1;
                elseif startsWith(dataType, "QF")
                    chainName = "PF-QF-" + extractAfter(dataType, "-");
                    dimIdx = 2;
                else
                    continue; % Unhandled type
                end
                
                % Find the corresponding chain object
                chainIdx = find([sensorChains.dataType] == chainName);
                if isempty(chainIdx)
                    continue; 
                end
                chainObj = sensorChains(chainIdx);
                
                % Extract GMM parameters
                mu     = chainObj.Means;                
                prob   = chainObj.prob;
                sigma2 = chainObj.Variances;
                
                % 2. Identify Clean and Noise Components
                [~, cleanIdx] = max(prob); % The healthy historical cluster
                [~, noiseIdx] = min(prob); % The cluster representing the corruption
                
                % Extract parameters for the specific dimension
               if chainName == "VM"
                    % Voltages are stable; use the tight clean variance
                    estData.Measurements(rowIdx) = mu(cleanIdx, 1);
                    estData.Variances(rowIdx)    = sigma2(1, 1, cleanIdx);
                else
                    % Power is dynamic; we must create the Phantom Anchor
                    
                    % 1. The Value: Safe Historical Mean (converted to p.u.)
                    zPseudo = mu(cleanIdx, dimIdx) / basVal;
                    if dataType == "PD" || dataType == "QD"
                        zPseudo = -abs(zPseudo); 
                    end
                    estData.Measurements(rowIdx) = zPseudo;
                    
                    % 2. The Variance: The massive GMM Noise Variance (converted to p.u.)
                    % This mathematically mimics the Huber down-weighting
                    varPseudo = sigma2(dimIdx, dimIdx, noiseIdx) / (basVal^2);
                    estData.Variances(rowIdx) = varPseudo; 
               end

                % 7. CRITICAL: Un-flag the row so the WLS solver keeps it in the Jacobian!
                % estData.Corrupt(rowIdx) = false;
            end
        end

        function estData = estCorruptData11(obj, data)
            estData = data;
            basVal = obj.Grid.baseVal; 
            
            % 1. Get the actual ROW INDICES of the corrupted measurements
            corruptRows = find(data.Corrupt);
            
            for k = 1:length(corruptRows)
                rowIdx = corruptRows(k);
                
                busID = data.Bus(rowIdx);
                dataType = string(data.DataType(rowIdx)); 
                
                % (Optional Safety) Skip PG/QG if prepareSensorMeas already merged them into PD/QD
                if dataType == "PG" || dataType == "QG"
                    continue;
                end
                
                % Get the sensor chains for this specific bus
                sensorChains = obj.sensorData(busID).Chains;
                if isempty(sensorChains)
                    continue;
                end
                
                % Map the DataType to the correct Chain Name and Dimension (1=P/V, 2=Q)
                if dataType == "VM"
                    chainName = "VM";
                    dimIdx = 1;
                elseif dataType == "PD"
                    chainName = "PD-QD";
                    dimIdx = 1;
                elseif dataType == "QD"
                    chainName = "PD-QD";
                    dimIdx = 2;
                elseif startsWith(dataType, "PF")
                    chainName = "PF-QF-" + extractAfter(dataType, "-");
                    dimIdx = 1;
                elseif startsWith(dataType, "QF")
                    chainName = "PF-QF-" + extractAfter(dataType, "-");
                    dimIdx = 2;
                else
                    continue; % Unhandled type
                end
                
                % Find the corresponding chain object
                chainIdx = find([sensorChains.dataType] == chainName);
                if isempty(chainIdx)
                    continue; 
                end
                chainObj = sensorChains(chainIdx);
                
                % Extract GMM parameters
                mu     = chainObj.Means;                
                prob   = chainObj.prob;
                sigma2 = chainObj.Variances;
                
                % 2. Find the "Clean" cluster (Highest Probability)
                [~, cleanIdx] = max(prob);
                
                % Extract Synthesized Value and Variance for the specific dimension
                if chainName == "VM"
                    zPseudo = mu(cleanIdx, 1);
                    varPseudo = sigma2(1, 1, cleanIdx);
                else
                    zPseudo = mu(cleanIdx, dimIdx);
                    varPseudo = sigma2(dimIdx, dimIdx, cleanIdx);
                end
                
                % 3. Apply Physical Scaling & WLS Conventions
                if chainName == "VM"
                    % VM is usually natively in p.u.
                    estData.Measurements(rowIdx) = zPseudo;
                    estData.Variances(rowIdx) = max(varPseudo, 1e-6); % Floor variance
                else
                    % Convert Power to per-unit
                    zPseudo = zPseudo / basVal;
                    varPseudo = varPseudo / (basVal^2);
                    
                    % MATPOWER WLS expects Loads (PD/QD) to be Negative Injections
                    if dataType == "PD" || dataType == "QD"
                        zPseudo = -abs(zPseudo); 
                    end
                    
                    estData.Measurements(rowIdx) = zPseudo;
                    
                    % Prevent Zero-Variance Matrix Explosion
                    estData.Variances(rowIdx) = max(varPseudo, 1e-6); 
                end
                
                % 4. CRITICAL: Un-flag the row so WLS actually uses the synthesized data!
                % estData.Corrupt(rowIdx) = false;
            end
        end

        function estData = estCorruptDataWrong(obj, data)
        
            % Initialize
            estData = data;
            basVal = obj.Grid.baseVal; % Get Base MVA for scaling
            
            % Get corrupted data 
            corruptBus = data.Bus(data.Corrupt);
            corruptBusDataType  = data.DataType(data.Corrupt);
            
            % Check if VM is there
            vmInd = double(any(ismember(corruptBusDataType, "VM")));
            
            % Compute the number of chains
            nChains = (numel(corruptBusDataType)+ vmInd)/2 ;
            
            % Go through all the measurements
            measCnt = 0;
            for iChain = 1:nChains
               
                % Map it to sensor processing chain
                sensorChains = obj.sensorData(corruptBus(iChain+measCnt)).Chains;
                chainsDataTypes = [sensorChains.dataType];
               
                chainIdx = contains(chainsDataTypes, corruptBusDataType(iChain+measCnt));
                corruptSensorChain = sensorChains(chainIdx);
                
                % Extract GMM parameters
                mu   = corruptSensorChain.Means;                
                prob = corruptSensorChain.prob;
                sigma2 = corruptSensorChain.Variances;
                
                % --- NEW: ISOLATE THE CLEAN COMPONENT ---
                % Find the component with the highest probability (The "Clean" cluster)
                [~, cleanIdx] = max(prob);
                
                % Extract strictly the clean mean and variance
                zPseudo_clean = mu(cleanIdx, :);
                
                if strcmp(corruptBusDataType(iChain+measCnt), "VM")
                    % 1D Data (Voltage)
                    % The variance is just the 1x1 scalar from the 3D Sigma matrix
                    var_clean = sigma2(1, 1, cleanIdx);
                    
                    % Modify corrupted data in the chain (VM is usually already p.u.)
                    estData.Measurements(corruptBus(iChain)) = zPseudo_clean;                    
                    estData.Variances(corruptBus(iChain)) = var_clean;
                    
                else
                    % 2D Data (P and Q)
                    % The variance is the diagonal of the 2x2 covariance matrix
                    var_clean = diag(sigma2(:,:,cleanIdx)); % Returns a 2x1 vector
                    
                    % Modify corrupted data in the chain
                    % CRITICAL FIX: Scale by basVal to convert to per-unit!
                    estData.Measurements(corruptBus(iChain:iChain+1)) = zPseudo_clean' / basVal;
                    
                    % CRITICAL FIX: Scale variance by basVal^2 !
                    estData.Variances(corruptBus(iChain:iChain+1)) = var_clean / (basVal^2);
                end
                
                % Measurement counter
                measCnt = measCnt + 1;
            end
        end
        function estData = estCorruptDataOld(obj, data)
        
            % Initialize
            estData = data;

            % Get corrupted data 
            corruptBus = data.Bus(data.Corrupt);
            corruptBusDataType  = data.DataType(data.Corrupt);

            % Check if VM is there
            vmInd = double(any(ismember(corruptBusDataType, "VM")));

            % Compute the number of chains
            nChains = (numel(corruptBusDataType)+ vmInd)/2 ;

            % Go through all the measurements
            measCnt = 0;
            for iChain = 1:nChains
               
                % Map it to sensor processing chain
                sensorChains = obj.sensorData(corruptBus(iChain+measCnt)).Chains;
                chainsDataTypes = [sensorChains.dataType];
               
                chainIdx = contains(chainsDataTypes, corruptBusDataType(iChain+measCnt));
                corruptSensorChain = sensorChains(chainIdx);

                % Exract mean and post prob.
                mu   = corruptSensorChain.Means;                
                prob = corruptSensorChain.prob;

                % Compute estimated data
                zPseudo = prob * mu;
                              
                if strcmp(corruptBusDataType(iChain+measCnt), "VM")
                    % Modify corrupted data in the chain
                    estData.Measurements(corruptBus(iChain)) = zPseudo';                    

                    % Compute corresponding variance
                    sigma2 = squeeze(corruptSensorChain.Variances);
                    estData.Variances(corruptBus(iChain)) = (prob * (mu - zPseudo).^2 + prob * sigma2)';

                else
                    % Modify corrupted data in the chain
                    estData.Measurements(corruptBus(iChain:iChain+1)) = zPseudo';

                    % Compute corresponding mean
                    sigma2 = corruptSensorChain.Variances;
                    sigma2PQcomp1 = diag(sigma2(:,:,1));
                    sigma2Qcomp2 = diag(sigma2(:,:,2));
                    sigma2PQ = [sigma2PQcomp1' ; sigma2Qcomp2'];

                    estData.Variances(corruptBus(iChain:iChain+1)) = (prob * (mu - zPseudo).^2 + prob * sigma2PQ)';
                end
                % Measurement counter
                measCnt = measCnt + 1;
            end
        end


        function W_iter = applyRobustWeights(~, W_base, raw_res, Options)
            % APPLYROBUSTWEIGHTS Adjusts the weight matrix for IRWLS estimators
            
            % Check if robust estimation is requested
            if isfield(Options, 'RobustEstimator') && strcmpi(Options.RobustEstimator, 'Huber')
                % Huber tuning constant (c = 1.5 is standard)
                huber_c = 1.5; 
                
                % 1. Compute Robust Scale (MAD)
                scale_s = 1.4826 * median(abs(raw_res)) + 1e-6;
                
                % 2. Standardize Residuals
                r_norm = raw_res / scale_s;
                
                % 3. Calculate Huber Weights
                huber_w = ones(size(r_norm)); % Default weight is 1
                outlier_idx = abs(r_norm) > huber_c;
                
                % Down-weight the outliers
                huber_w(outlier_idx) = huber_c ./ abs(r_norm(outlier_idx)); 
                
                % 4. Apply weights to the diagonal of the Base Weight Matrix
                W_iter = W_base * diag(huber_w);
            else
                % If ordinary WLS, just return the base matrix unchanged
                W_iter = W_base;
            end
        end

        function updateMCProgress(~, data, nSim)
            % If the worker sent a string (e.g., a warning), print it
            if ischar(data)
                fprintf('%s\n', data);
            % If the worker sent an iteration number, calculate progress
            elseif isnumeric(data)
                % Print a status update every 5% to avoid console spam
                updateInterval = max(1, floor(nSim / 20)); 
                if mod(data, updateInterval) == 0
                    fprintf('Completed %d / %d simulations (%.0f%%)\n', ...
                        data, nSim, (data/nSim)*100);
                end
            end
        end

        function W = computeWeightMatrix(~, zTable, commVarTable, Options)
            
            zVar = zTable.Variances;
            commsVar = commVarTable.Comm_Variance;
            if strcmp(Options.weightMatrix, 'unity')
                W = diag(ones(length(zVar),1));

            elseif strcmp(Options.weightMatrix, 'comms')
                W = diag(1./commsVar);

            elseif strcmp(Options.weightMatrix, 'meas-comms')
                W = diag(1./(zTable.Variances + commsVar) );
            else
                error('Unknown weight matrix type!')
            end

            if ~(Options.FilteringOverride)
                % Remove rows corresponding to corrupted data
                excludeFlag = xor(zTable.Corrupt, zTable.Override);
                excludeFlag = logical(excludeFlag );
                W(excludeFlag,:) = [];
                W(:, excludeFlag) = [];
            end

        end
        
        function prepareDecoupledMat(obj)
            % Prepare the Bp and Bpp required for the decoupled power
            % initialisation algorithm.
            % Prapering Bp matrix: Remove slack bus to from the admittance
            % matrix.
            % Preparing Bpp matrix: Remove the slack and PV buses from the
            % admittance matrix.

            % Identify the slack and PV buses
            busType = {obj.Grid.Bus.Type};

            slackBusInd = find(string(busType) == 'Ref');
            pvBusInd = find(string(busType) == 'PV');

            % Remove slack bus row and column from Bp
            Bp = obj.Grid.Bp;
            Bp(slackBusInd,:) = [];
            Bp(:,slackBusInd) = [];
            obj.Grid.Bp = full(Bp);

            % Remove the slack and PV buses
            Bpp = obj.Grid.Bpp;
            Bpp(slackBusInd,:) = [];
            Bpp(:,slackBusInd) = [];
            Bpp(pvBusInd,:) = [];
            Bpp(:,pvBusInd) = [];
            obj.Grid.Bpp = full(Bpp);


        end

        function zTable = prepareSensorMeas(obj, busMapping, Idx)
            % Inline function to split data type into its elements

            % Grid base value
            basVal = obj.Grid.baseVal;
            % Idx = 1; % data index
            nBus = numel(obj.Grid.Bus);

            % Get all data types in chain connected to iBus bus
            chainDataTypeCell = cell(1,nBus);
            for iBus = 1:nBus
                nChains = numel(obj.sensorData(iBus).Chains);
                chainDataType = [];
                for iChains = 1:nChains
                    chainDataType = [chainDataType, ...
                        obj.sensorData(iBus).Chains(iChains).dataType];
                end
                chainDataTypeCell{iBus} = chainDataType;
            end

            % Initialize table
            zTable = busMapping;

            % Go through all the measurements/sensors
            nMeas = height(zTable);

            for iMeas =1:nMeas
                % Look bus number
                iBus = zTable.Bus(iMeas);

                % Look up measurement data type
                datatype = zTable.DataType(iMeas);

                % datatype = SensorConnector.DataType(dataTypeInd);
                branchStr = "-"+ string(zTable.Branch(iMeas));
                branchStr(ismissing(branchStr)) = "";

                dataType = datatype + branchStr;

                % Find corresponding chain
                dataTypeChainIdx = contains(chainDataTypeCell{iBus}, datatype) & ...
                    contains(chainDataTypeCell{iBus}, branchStr);
                selectedChain =  obj.sensorData(iBus).Chains(dataTypeChainIdx);

                if strcmp(dataType, "VM")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,:);
                    zTable.Variances(iMeas) = selectedChain.Variances(:,:,1);

                elseif strcmp(dataType, "PG") || strcmp(dataType, "PD") || strcmp(extractBetween(dataType,1,2), "PF")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,1);
                    zTable.Variances(iMeas) = selectedChain.Variances(1,1,1);

                elseif strcmp(dataType, "QG") || strcmp(dataType, "QD") || strcmp(extractBetween(dataType,1,2), "QF")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,2);
                    zTable.Variances(iMeas) = selectedChain.Variances(2,2,1);
                end

                zTable.Corrupt(iMeas) = selectedChain.Corrupt(Idx);
            end

            % ---------------------------------------------------------
            % STEP 1: GLOBAL NEGATION FOR LOADS
            % Solver expects Net Injection (Gen - Load).
            % Sensors report Load as Positive. We must flip sign to Negative.
            % ---------------------------------------------------------
            pdRows = (zTable.DataType == "PD");
            zTable.Measurements(pdRows) = -zTable.Measurements(pdRows);

            qdRows = (zTable.DataType == "QD");
            zTable.Measurements(qdRows) = -zTable.Measurements(qdRows);

            % ---------------------------------------------------------
            % STEP 2: MERGE GENERATION (PG) INTO NET INJECTION (PD)
            % ---------------------------------------------------------
            pgRowsInd = find(zTable.DataType == "PG");
            for k = 1:numel(pgRowsInd)
                pgIdx = pgRowsInd(k);
                busID = zTable.Bus(pgIdx);

                % Check if this bus also has a PD row
                pdIdx = find(zTable.Bus == busID & zTable.DataType == "PD");

                if ~isempty(pdIdx)
                    % Add Generation (PG) to the already negated Load (-PD)
                    % Result: Net Injection = PG - PD
                    zTable.Measurements(pdIdx) = zTable.Measurements(pdIdx) + zTable.Measurements(pgIdx);

                    % MARK PG AS CORRUPT so WLS ignores it (prevents conflict)
                    zTable.Corrupt(pgIdx) = true;
                end
            end

            % Repeat for Reactive Power (QG -> QD)
            qgRowsInd = find(zTable.DataType == "QG");
            for k = 1:numel(qgRowsInd)
                qgIdx = qgRowsInd(k);
                busID = zTable.Bus(qgIdx);

                qdIdx = find(zTable.Bus == busID & zTable.DataType == "QD");

                if ~isempty(qdIdx)
                    zTable.Measurements(qdIdx) = zTable.Measurements(qdIdx) + zTable.Measurements(qgIdx);
                    zTable.Corrupt(qgIdx) = true;
                end
            end

            % ---------------------------------------------------------
            % STEP 3: NORMALIZE TO PER-UNIT
            % ---------------------------------------------------------
            nMeas = height(zTable);
            for iMeas = 1:nMeas
                % Only divide Power quantities, not Voltage Angles or VM (if VM is already p.u.)
                if ~strcmp(zTable.DataType(iMeas), "VM") && ~strcmp(zTable.DataType(iMeas), "VA")
                    zTable.Measurements(iMeas) = zTable.Measurements(iMeas) / basVal;
                end
            end

            % Add override column initialized to false
            zTable.Override = false(nMeas,1);

        end

        function zTable = prepareSensorMeasOld(obj, busMapping)

            % Inline function to split data type into its elements
            % splitPair = @(s) split(s, "-")';

            % Grid base value
            basVal = obj.Grid.baseVal;

            % !! Data index set to 1 for testing. Account for MC
            % simulation
            Idx = 1; % data index
            nBus = numel(obj.Grid.Bus);

            % Get all data types in chain connected to iBus bus
            chainDataTypeCell = cell(1,nBus);

            for iBus = 1:nBus
                nChains = numel(obj.sensorData(iBus).Chains);
                chainDataType = [];
                for iChains = 1:nChains
                    chainDataType = [chainDataType, ...
                        obj.sensorData(iBus).Chains(iChains).dataType];
                end
                chainDataTypeCell{iBus} = chainDataType;
            end


            % Initialize cell array
            % z = {};
            % Initialize table
            zTable = busMapping;

            % Go through all the measurements/sensors
            nMeas = height(zTable);

            % Table counter
            iMeas = 1;
            while (iMeas <=nMeas)
                % Look bus number
                iBus = zTable.Bus(iMeas);

                % Look up measurement data type
                dataType = zTable.DataType(iMeas);

                % Find corresponding chain
                dataTypeChainIdx = contains(chainDataTypeCell{iBus}, dataType);
                selectedChain =  obj.sensorData(iBus).Chains(dataTypeChainIdx);

                % Using above bus # and data type lookup the rest of values
                % such as measurement, variance and corruption.
                % Data is stored in the chains as:
                % VM - 1D.
                % PG-QG, PD-QD, PF-QF: 2D.

                if strcmp(dataType, "VM")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,:);
                    zTable.Variances(iMeas) = selectedChain.Variances(1);

                elseif strcmp(dataType, "PG") || strcmp(dataType, "PD") || strcmp(dataType, "PF")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,1);
                    zTable.Variances(iMeas) = selectedChain.Variances(1);

                elseif strcmp(dataType, "QG") || strcmp(dataType, "QD") || strcmp(dataType, "QF")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,2);
                    zTable.Variances(iMeas) = selectedChain.Variances(2);

                end

                zTable.Corrupt(iMeas) = selectedChain.Corrupt(Idx);

                % Increment measurement counter
                iMeas = iMeas + 1;
            end


            % nBus = numel(obj.sensorData);
            % for iBus = 1:nBus
            %     nChains = numel(obj.sensorData(iBus).Chains);
            %     for iChains = 1:nChains
            %
            %         % Account for VM data being 1D
            %         if strcmp(obj.sensorData(iBus).Chains(iChains).dataType ,"VM")
            %             r = 1;
            %         else
            %             r = 2;
            %         end
            %
            %         % Collect data in a cell array
            %         z(end + 1, :) = {obj.sensorData(iBus).Bus_ID*ones(r,1), ...
            %             splitPair(obj.sensorData(iBus).Chains(iChains).dataType), ...
            %             obj.sensorData(iBus).Chains(iChains).data(Idx,:)', ...
            %             obj.sensorData(iBus).Chains(iChains).Variances, ...
            %             obj.sensorData(iBus).Chains(iChains).Corrupt(Idx)*ones(r,1)};
            %     end
            % end
            % % Convert cell array to table
            % zTable = table(vertcat(z{:,1}), ...
            %     vertcat(z{:,2}), ...
            %     vertcat(z{:,3}), ...
            %     vertcat(z{:,4}), ...
            %     vertcat(z{:,5}), ...
            %     'VariableNames',  ...
            %     {'Bus','DataType','Measurements','Variances', 'Corrupt'});

            % Combine PG,QG with PD,QD
            % Find the PG/QG elements
            % pgInx = (zTable{:,2} == "PG");
            % qgInx = (zTable{:,2} == "QG");
            pgInx = (zTable.DataType(:) == "PG");
            qgInx = (zTable.DataType(:) == "QG");

            % Get corresponding buses
            % pgBus = zTable{pgInx,1};
            % qgBus = zTable{qgInx,1};
            pgBus = zTable.Bus(pgInx,:);
            qgBus = zTable.Bus(qgInx,:);

            % Check if there are PD
            for iBus = 1:numel(pgBus)
                % busIdx = (zTable{:,1} == pgBus(iBus));
                busIdx = (zTable.Bus == pgBus(iBus));

                if any(zTable.DataType(busIdx,:) == "PD")
                    % Combine if there are (PG - PD)
                    % pdIdx = logical(busIdx.*(zTable{:,2} == "PD"));
                    % pgIdx = logical(busIdx.*(zTable{:,2} == "PG"));
                    % zTable{pdIdx,3} = (zTable{pgIdx,3} - zTable{pdIdx,3});
                    pdIdx = logical(busIdx.*(zTable.DataType(:) == "PD"));
                    pgIdx = logical(busIdx.*(zTable.DataType(:) == "PG"));
                    zTable{pdIdx,3} = (zTable.Measurements(pgIdx,:) - zTable.Measurements(pdIdx,:));

                    % FIX: Mark the original PG row as Corrupt so it is ignored by WLS
                    zTable{pgIdx, 5} = 1;
                end
            end

            % Check if there are QD
            for iBus = 1:numel(qgBus)
                % busIdx = (zTable{:,1} == qgBus(iBus));
                busIdx = (zTable.Bus == qgBus(iBus));

                if any(zTable.DataType(busIdx,:) == "QD")
                    % Combine if there are (QG - QD)
                    qdIdx = logical(busIdx.*(zTable.DataType(:) == "QD"));
                    qgIdx = logical(busIdx.*(zTable.DataType(:) == "QG"));
                    zTable{qdIdx,3} = (zTable.Measurements(qgIdx,:) - zTable.Measurements(qdIdx,:));

                    % FIX: Mark the original QG row as Corrupt
                    zTable{qgIdx, 5} = 1;
                end
            end

            % Normalize power values
            nMeas = height(zTable);

            for iMeas = 1:nMeas
                if ~strcmp(zTable.DataType(iMeas), "VM")
                    zTable.Measurements(iMeas) = zTable.Measurements(iMeas) / basVal;
                end
            end

        end

        function busTypeInd = getBusTypeInd(obj, busType)

            buses = {obj.Grid.Bus.Type};
            busTypeInd = find(string(buses) == busType);

        end

        function [Pinj, Qinj] = collectPQinj(obj, data)

            slackBusInd = obj.getBusTypeInd('Ref');
            pvBusInd = obj.getBusTypeInd('PV');

            Ptot = obj.extractMeas(data,"PD");
            Qtot = obj.extractMeas(data,"QD");

            % Replace unavailable values with nominal values
            % Pnominal = abs(obj.Grid.TrueNetPower.PD);
            % Qnominal = abs(obj.Grid.TrueNetPower.QD);
            Pnominal = obj.Grid.TrueNetPower.PD/ obj.Grid.baseVal;
            Qnominal = obj.Grid.TrueNetPower.QD/ obj.Grid.baseVal;

            pNanId = isnan(Ptot);
            Ptot(pNanId) = Pnominal(pNanId);

            qNanId = isnan(Qtot);
            Qtot(qNanId) = Qnominal(qNanId);

            % --- FIX: REMOVE "/ obj.Grid.baseVal" ---
            % The input 'data' (zTable) is already in per-unit.
            Pinj = Ptot;
            Pinj(slackBusInd) = [];

            Qinj = Qtot;
            Qinj(slackBusInd) = [];
            Qinj(pvBusInd) = [];

        end

        function [Pinj, Qinj] = collectPQinjOld(obj, data)


            slackBusInd = obj.getBusTypeInd('Ref');
            pvBusInd = obj.getBusTypeInd('PV');

            Ptot = obj.extractMeas(data,"PD");
            Qtot = obj.extractMeas(data,"QD");



            % Replace unavailable values with nominal values
            Pnominal = abs(obj.Grid.TrueNetPower.PD);
            Qnominal = abs(obj.Grid.TrueNetPower.QD);

            pNanId = isnan(Ptot);
            Ptot(pNanId) = Pnominal(pNanId);

            qNanId = isnan(Qtot);
            Qtot(pNanId) = Qnominal(qNanId);

            % Pinj = Ptot / obj.Grid.baseVal;
            Pinj = Ptot;
            % Remove slack bus value
            Pinj(slackBusInd) = [];

            % Qinj = Qtot / obj.Grid.baseVal;
            Qinj = Qtot;
            % Remove slack and OV buses values
            Qinj(slackBusInd) = [];
            Qinj(pvBusInd) = [];



        end

        function commVarTable = prepareCommVar(obj, Idx)
            %
            % Inline function to split data type into its elements
            splitPair = @(s) split(s, "-");

            commVarCell = {};
            % Idx = 1; % data index
            nBus = numel(obj.sensorData);
            for iBus = 1:nBus
                nChains = numel(obj.sensorData(iBus).Chains);
                for iChains = 1:nChains

                    if strcmp(obj.sensorData(iBus).Chains(iChains).dataType ,"VM")
                        r = 1;
                    else
                        r = 2;
                    end

                    commVarCell(end + 1, :) = {obj.sensorData(iBus).Bus_ID*ones(r,1), ...
                        splitPair(obj.sensorData(iBus).Chains(iChains).dataType), ...
                        obj.sensorData(iBus).Chains(iChains).commVar(Idx,:)'};
                    
                    % Account for the flow measurement
                    if size(commVarCell{end, 2},1) == 3
                        intr = commVarCell{end, 2};
                        intr(1) = intr(1) + "-" + intr(3);
                        intr(2) = intr(2) + "-" + intr(3);
                        intr(3) = [];
                        commVarCell{end, 2} = intr;
                    end
                end
            end
            % Convert to table
            commVarTable = table(vertcat(commVarCell{:,1}), ...
                vertcat(commVarCell{:,2}), ...
                vertcat(commVarCell{:,3}), ...
                'VariableNames', {'Bus','DataType','Comm_Variance'});

            % Combine PG,QG with PD,QD
            % Find the PG/QG elements
            pgInx = (commVarTable{:,2} == "PG");
            qgInx = (commVarTable{:,2} == "QG");

            % Get corresponding buses
            pgBus = commVarTable{pgInx,1};
            qgBus = commVarTable{qgInx,1};

            % Check if there are PD
            for iBus = 1:numel(pgBus)
                busIdx = (commVarTable{:,1} == pgBus(iBus));

                if any(commVarTable{busIdx,2} == "PD")
                    % Combine if there are (PG - PD)
                    pdIdx = logical(busIdx.*(commVarTable{:,2} == "PD"));
                    pgIdx = logical(busIdx.*(commVarTable{:,2} == "PG"));
                    commVarTable{pdIdx,3} = commVarTable{pgIdx,3} + commVarTable{pdIdx,3};

                end
            end

            % Check if there are QD
            for iBus = 1:numel(qgBus)
                busIdx = (commVarTable{:,1} == qgBus(iBus));

                if any(commVarTable{busIdx,2} == "QD")
                    % Combine if there are (QG - QD)
                    qdIdx = logical(busIdx.*(commVarTable{:,2} == "QD"));
                    qgIdx = logical(busIdx.*(commVarTable{:,2} == "QG"));
                    commVarTable{qdIdx,3} = commVarTable{qgIdx,3} + commVarTable{qdIdx,3};

                end
            end


        end

        function zHatTable = computeMeas(obj, x, z)

            define_constants;
            mpc  = obj.Grid.mpc;
            Ybus = obj.Grid.Ybus;
            Yf   = obj.Grid.Yf;

            % Reconstruct Complex Voltage Vector
            V = x.VM .* exp(1i * deg2rad(x.VA) );

            % Injection power values

            % Compute Power Injections (S = P + jQ)
            S_inj = V .* conj(Ybus * V);

            % Compute Power Flows
            S_flow_f = V(mpc.branch(:, 1)) .* conj(Yf * V);

            % Get Real and Reactive Powers
            P_inj_all = real(S_inj);
            Q_inj_all = imag(S_inj);

            P_flow_f_all = real(S_flow_f);
            Q_flow_f_all = imag(S_flow_f);


            % Select the Specific Measurements
            nData = height(z);
            zHat  = zeros(nData,1);

            for iData = 1:nData
                idx = double(z{iData, 'Bus'});
                switch char(z{iData, 'DataType'})

                    case 'VM'
                        % Voltage Magnitude (Pass-through from State)
                        zHat(iData) = x.VM(idx);

                    case 'VA'
                        % Voltage Angle (Pass-through from State)
                        % Note: Returns Degrees because x.VA is in Degrees
                        zHat(iData) = x.VA(idx);
                    case {'PD', 'PG'}
                        zHat(iData) = P_inj_all(idx);

                    case {'QD','QG'}
                        zHat(iData) = Q_inj_all(idx);

                    case 'PF'
                        % zHat(iData) = P_flow_f_all(idx);
                        branchIdx = double(z{iData, 'Branch'});
                        zHat(iData) = P_flow_f_all(branchIdx);

                    case 'QF'
                        % zHat(iData) = Q_flow_f_all(idx);
                        branchIdx = double(z{iData, 'Branch'});
                        zHat(iData) = Q_flow_f_all(branchIdx);

                end

            end

            zHatTable = z;
            zHatTable.Measurements = zHat;

            % Add Override column
            zHatTable.Override = false(nData, 1);

            % zHatTable = addvars(zHatTable, zHat, 'NewVariableNames', 'Measurements');

        end

        function r = computeResiduals(obj, zTable, zHatTable)

            nMeas = height(zTable);

            % Precallocate residual table
            r = zTable(:,1:2);
            r.Residuals = zeros(height(r),1);

            for iMeas = 1:nMeas
                if ~strcmp(table2array(zTable(iMeas,2)), table2array(zHatTable(iMeas,2)))
                    error('Measurements mismatch!')
                else
                    r.Residuals(iMeas) = zTable.Measurements(iMeas) ...
                        - zHatTable.Measurements(iMeas);
                end
            end

            % Remove corrupted data
            if all(zTable.Override)
                excludeFlag = false(nMeas,1);
            else
                excludeFlag = xor(zTable.Corrupt, zTable.Override);
                r(excludeFlag,:) = [];
            end

        end

        function H = computeJacobian(obj, x, measMapping)
            % COMPUTEJACOBIAN Computes H matrix for WLS State Estimation
            % UPDATED: Uses 'Branch' column for flow sensors and adds PT/QT support.

            % 1. Setup
            mpc = obj.Grid.mpc;
            Yf = obj.Grid.Yf;
            Yt = obj.Grid.Yt;
            Ybus = obj.Grid.Ybus;

            % 2. Reconstruct Complex Voltage Vector
            % Note: x.VA is in degrees, convert to radians
            V = x.VM .* exp(1i * deg2rad(x.VA));

            % 3. Compute Partial Derivatives (The "Building Blocks")
            % A) INJECTIONS (dSbus_dV)
            % This returns the complex derivatives dS/dVa and dS/dVm (Size: nBus x nBus)
            [dSbus_dVa, dSbus_dVm] = dSbus_dV(Ybus, V);

            % B) FLOWS (dSbr_dV)
            % This returns complex derivatives for From (f) and To (t) ends (Size: nBranch x nBus)
            [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, ~, ~] = dSbr_dV(mpc.branch, Yf, Yt, V);

            % 4. Assemble H Matrix Row-by-Row
            nMeas = height(measMapping);
            nBus  = size(mpc.bus, 1);

            % Preallocate H (Sparse is faster for large grids)
            % Columns: [Angle_1 ... Angle_N, Voltage_1 ... Voltage_N]
            H = sparse(nMeas, 2*nBus);

            % corruptVec = measMapping.Corrupt;
            if all(measMapping.Override)
                corruptVec = false(height(measMapping),1);
            else
                corruptVec = xor(measMapping.Corrupt, measMapping.Override);
            end

            for i = 1:nMeas
                type = char(measMapping.DataType(i));

                % --- FIX: SELECT CORRECT INDEX SOURCE ---
                % For Flows: Use the Branch ID (e.g., Branch 15)
                % For Nodes: Use the Bus ID (e.g., Bus 8)
                if ismember(type, {'PF', 'QF', 'PT', 'QT'})
                    idx = measMapping.Branch(i);
                else
                    idx = measMapping.Bus(i);
                end

                corruptFlag = corruptVec(i);

                % If the measurement is corrupted leave the row zero and move on
                if corruptFlag
                    continue
                end

                switch type
                    case {'PD', 'PG'} % Active Injection (P)
                        H(i, 1:nBus)     = real(dSbus_dVa(idx, :)); % dP/dVa
                        H(i, nBus+1:end) = real(dSbus_dVm(idx, :)); % dP/dVm

                    case {'QD', 'QG'} % Reactive Injection (Q)
                        H(i, 1:nBus)     = imag(dSbus_dVa(idx, :)); % dQ/dVa
                        H(i, nBus+1:end) = imag(dSbus_dVm(idx, :)); % dQ/dVm

                    case 'PF' % Active Flow (From Bus -> To Bus)
                        H(i, 1:nBus)     = real(dSf_dVa(idx, :));
                        H(i, nBus+1:end) = real(dSf_dVm(idx, :));

                    case 'QF' % Reactive Flow (From Bus -> To Bus)
                        H(i, 1:nBus)     = imag(dSf_dVa(idx, :));
                        H(i, nBus+1:end) = imag(dSf_dVm(idx, :));

                        % --- NEW CASES FOR RADIAL BUSES ---
                    case 'PT' % Active Flow (To Bus -> From Bus)
                        % Uses dSt (To-side derivatives)
                        H(i, 1:nBus)     = real(dSt_dVa(idx, :));
                        H(i, nBus+1:end) = real(dSt_dVm(idx, :));

                    case 'QT' % Reactive Flow (To Bus -> From Bus)
                        % Uses dSt (To-side derivatives)
                        H(i, 1:nBus)     = imag(dSt_dVa(idx, :));
                        H(i, nBus+1:end) = imag(dSt_dVm(idx, :));

                    case 'VM' % Voltage Magnitude
                        % Derivative is 0 for angles, 1 for the specific bus magnitude
                        H(i, nBus+idx) = 1;

                    case 'VA' % Voltage Angle
                        % Derivative is 1 for the specific bus angle, 0 for magnitudes
                        H(i, idx) = 180/pi;
                end
            end

            % Remove corrupted rows
            H(corruptVec,:) = [];
        end

        function meas = extractMeas(obj, measurements, dataType)


            % Go through the measurements and extract PD
            nBus = numel(obj.Grid.Bus);
            meas = nan(nBus,1);
            measIdx = (measurements.DataType == dataType );
            busMeas = measurements.Bus(measIdx);
            corruptedBusInd = logical(measurements.Corrupt(busMeas));
            Val = measurements.Measurements(measIdx) ;
            meas(busMeas) = Val;
            meas(corruptedBusInd) = NaN;

        end

        function z = filterOverride(obj, x, FilteringOverride)

            z = x;
            if FilteringOverride
                z.Override = true(height(x),1);
            else
                nBus = numel(obj.Grid.Bus);
               
                % Check each bus for corrupted data
                for iBus = 1:nBus
                    corruptBusIdx = (x.Bus == iBus);
                    if all( x.Corrupt(corruptBusIdx) )                      
                        z.Override(corruptBusIdx) = true;
                        continue
                    end

                end
            end
        end

    end
end