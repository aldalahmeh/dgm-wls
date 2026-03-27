%%
classdef GridEdgeNode < handle
    % GRIDEDGENODE - Representation of a grid edge node used in simulation
    %   OBJ = GRIDEDGENODE() encapsulates simulation parameters, the grid
    %   structure, sensor data, and bookkeeping for completed sims.
    %
    %   Properties:
    %       nSim            - number of simulations to run
    %       Grid            - grid topology / data structure
    %       sensorData      - collected or input sensor measurements
    %       gridSensorNoiseStd - sensor noise standard deviation
    %       wlsParam        - parameters for weighted least squares
    %       simsCompleted   - count of completed simulations
    %       DEBUG           - flag to enable debug behavior
    % GRIDEDGENODE Summary of this class goes here
    %   Detailed explanation goes here
    properties
        nSim
        Grid
        sensorData
        gridSensorNoiseStd
        wlsParam
        simsCompleted
        DEBUG
    end

    methods
        %% Constructor
        function obj = GridEdgeNode(Grid, nSim, gridSensorNoiseStd, maxIter, tol)
          % GRIDEDGENODE - Construct GridEdgeNode object
          %
          % Input arguments:
          % Grid                - grid data structure containing Bus field
          % nSim                - number of simulations
          % gridSensorNoiseStd  - sensor noise standard deviation
          % maxIter (opt)       - WLS max iterations
          % tol (opt)           - WLS convergence tolerance
          %
          % Output arguments:
          % obj - constructed GridEdgeNode instance
            %GridEdgeNode Construct an instance of this class
            %   Detailed explanation goes here
            % Debug flag
            obj.DEBUG = false;
            obj.nSim = nSim;
            obj.Grid = Grid;
            obj.gridSensorNoiseStd = gridSensorNoiseStd;
            % WLS parameters
            if nargin == 5
                obj.wlsParam.maxIter = maxIter;
                obj.wlsParam.tol = tol;
            elseif  nargin == 4
                obj.wlsParam.maxIter = maxIter;
                obj.wlsParam.tol = 1e-4;
            elseif nargin < 4
                obj.wlsParam.maxIter = 20;
                obj.wlsParam.tol = 1e-4;
            end
            obj.wlsParam.alpha = 1;
            % Prepare the Bp matrix
            obj.prepareDecoupledMat();
            nBus = numel(Grid.Bus);
            % Preallocate sensorData struct
            obj.sensorData = repmat(struct('Bus_ID', [], 'Chains', []), 1, nBus );
        end
        %% Getters
        function datatype = getChainDataType(obj, BusNum, ChainNum)
          % GETCHAINDATATYPE - Return datatype for a specific chain on a bus
            datatype = obj.sensorData(BusNum).Chains(ChainNum).dataType;
        end
        function busID = getBusID(obj, iBus)
          % GETBUSID - Return the stored bus identifier for a given index
            busID = obj.sensorData(iBus).Bus_ID;
        end
        function gridSensorNoiseStd = getGridSensorNoiseStd(obj)
          % GETGRIDSENSORNOISESTD - Access global grid sensor noise standard dev
            gridSensorNoiseStd = obj.gridSensorNoiseStd;
        end
        function  commVar = getCommVar(obj, iBus, iChains)
          % GETCOMMVAR - Return communication variance for a specific chain
            commVar = obj.sensorData(iBus).Chains(iChains).commVar;
        end
        %%
        function receive(obj, rxData)
          % RECEIVE - Update stored sensorData entry with incoming packet
            % Get sensor ID
            sensorID = rxData.Bus_ID;
            obj.sensorData(sensorID) = rxData;
        end

        function initEst = computeInitEst(obj, initType, measurements )
          % COMPUTEINITEST - Generate initial voltage magnitude and angle estimates
          %
          % Input arguments:
          % obj          - estimator object containing Grid and helpers
          % initType     - 'flat' or 'warm' initialization strategy
          % measurements - measurement set used for 'warm' initialization
          %
          % Output arguments:
          % initEst      - table with Bus_ID, VM (voltage magnitude), VA (angle)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            nBus = numel(obj.Grid.Bus);
            if strcmpi(initType,'flat')
                vEstArray = ones(nBus,1);
                vAngEstArray = zeros(nBus, 1);
            elseif strcmpi(initType,'warm')
                % Collect measured active and reactive power injections
                [Pinj, Qinj] = obj.collectPQinj(measurements);
                % Get susceptance-like matrices for linearized relations
                Bp  = obj.Grid.getBp();
                Bpp = obj.Grid.getBpp();
                % Compute linearized initial voltage angles (degrees)
                vAngEst = rad2deg(-Bp \ Pinj);
                % Append a zero for the reference (slack) bus value
                nBus = numel(obj.Grid.Bus);
                slackBusInd = obj.getBusTypeInd('Ref');
                nonSlackIdx = setdiff(1:nBus, slackBusInd);
                vAngEstArray = zeros(nBus, 1);
                vAngEstArray(nonSlackIdx) = vAngEst;
                % Compute the linearized initial voltage magnitudes
                vEst  = 1.0 + Bpp \ Qinj;
                vEstArray = ones(nBus,1);
                pqBusInd = obj.getBusTypeInd('PQ');
                vEstArray(pqBusInd) = vEst;
            else
                error('Undefined state initialization type!')
            end
            % Pack results into a table for downstream processing
            initEst = table((1:nBus)', vEstArray, vAngEstArray,...
                'VariableNames',{'Bus_ID', 'VM', 'VA'});
        end

        function [stateEst, convergenceData] = computeStateEst(obj, Options)
          % COMPUTESTATEEST - Run Monte Carlo state estimation over nSim realizations
          %
          % Input arguments:
          % obj     - estimator object containing grid, parameters, and helpers
          % Options - struct controlling execution (e.g., DebugMode)
          %
          % Output arguments:
          % stateEst       - 1xnSim cell array of state estimates
          % convergenceData- 1xnSim vector of convergence iterations
            % Collect parameters
            busMeasMapping = obj.Grid.sensorTable;
            maxIter = obj.wlsParam.maxIter;
            tol     = obj.wlsParam.tol;
            alpha   = obj.wlsParam.alpha;
            nSim    = obj.nSim;
            % Setup DataQueue once (client-side)
            obj.simsCompleted = 0;
            q = parallel.pool.DataQueue;
            afterEach(q, @(data) obj.updateMCProgress(data, nSim));
            % Decide pool size and reuse existing pool if present
            useDebug = isfield(Options, 'DebugMode') && Options.DebugMode == true;
            if useDebug
                % Sequential execution: ensure we do not create a pool
                fprintf('--- RUNNING IN DEBUG MODE (Sequential) ---\n');
                numWorkers = 0;
            else
                % Parallel execution: reuse an existing pool or create one if absent
                desiredWorkers = 20;
                pool = gcp('nocreate');
                if isempty(pool)
                    % Create a pool only if none exists
                    pool = parpool(desiredWorkers);
                else
                    % If a pool exists but size differs, optionally warn or try to adjust
                    if pool.NumWorkers ~= desiredWorkers
                        fprintf('Existing pool with %d workers found (desired %d). Reusing existing pool.\n', ...
                            pool.NumWorkers, desiredWorkers);
                        % Optionally: you could shut and recreate here, but to keep pool alive we reuse it.
                    end
                end
                numWorkers = pool.NumWorkers;
                fprintf('Starting Monte Carlo simulation (%d runs) on %d workers...\n', nSim, numWorkers);
            end
            % Preallocate outputs
            stateEst = cell(1, nSim);
            convergenceData = zeros(1, nSim);
            % Choose loop type based on DebugMode
            if useDebug
                % Run serially for easier debugging and deterministic ordering
                for iSim = 1:nSim
                    [stateEst{iSim}, convergenceData(iSim)] = computeSingleStateEst(obj, iSim, ...
                        busMeasMapping, Options, q, maxIter, tol, alpha);
                end
            else
                % Parallel loop: each iteration is independent Monte Carlo run
                parfor iSim = 1:nSim
                    [sEst, convIter] = computeSingleStateEst(obj, iSim, ...
                        busMeasMapping, Options, q, maxIter, tol, alpha);
                    stateEst{iSim} = sEst;
                    convergenceData(iSim) = convIter;
                end
            end
        end

        function [stateOut, convIter] = computeSingleStateEst(obj, iSim, busMeasMapping, Options, q, maxIter, tol, alpha)
          % COMPUTESINGLESTATEEST - Perform one WLS state estimation Monte-Carlo run
          %
          % Input arguments:
          % obj - estimator object
          % iSim - simulation/run index
          % busMeasMapping - mapping from buses to measurements
          % Options - struct of algorithm options
          % q - DataQueue for progress/errors
          % maxIter - maximum WLS iterations
          % tol - convergence tolerance on state update
          % alpha - damping/step-length factor
          %
          % Output arguments:
          % stateOut - estimated bus state (best effort if not converged)
          % convIter - iteration at which convergence occurred (0 if none)
            % One Monte-Carlo simulation run. Must be a subfunction (not nested) so parfor can call it.
            wlsConverged = false;
            stateOut = [];       % best effort result
            convIter = 0;
            % Prepare measurements
            zTable = obj.prepareSensorMeas(busMeasMapping, iSim);
            zTable = obj.filterOverride(zTable, Options.FilteringOverride);
            % Initial estimate
            if strcmp(Options.WlsInit, 'warm')
                xCurrent = obj.computeInitEst('Warm', zTable);
            elseif strcmp(Options.WlsInit, 'flat')
                xCurrent = obj.computeInitEst('Flat', zTable);
            else
                error('Unknown initialization method!')
            end
            % Base weight
            W_base = obj.computeWeightMatrix(zTable, Options);
            nBus = height(xCurrent);
            for iter = 1:maxIter
                % Compute predicted measurements and residuals
                zHatTable = obj.computeMeas(xCurrent, busMeasMapping);
                zHatTable.Corrupt = zTable.Corrupt;
                zHatTable.Override = zTable.Override;
                % Estimate corrupted data if enabled
                if ~strcmp(Options.dataImpute,'None') && any(zTable.Corrupt)
                    zTable = obj.estCorruptData(zTable, Options);
                end
                r = obj.computeResiduals(zTable, zHatTable);
                % Robust weights (or returns W_base)
                W_iter = obj.applyRobustWeights(W_base, r.Residuals, Options);
                % Jacobian and reductions
                H = obj.computeJacobian(xCurrent, zHatTable);
                Hred = H;
                Hred(:,1) = [];
                G = Hred' * W_iter * Hred;
                rhs = Hred' * W_iter * r.Residuals;
                try
                    lastwarn('');
                    dx_red = G \ rhs;
                    [~, msgid] = lastwarn;
                    if strcmp(msgid, 'MATLAB:nearlySingularMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
                        error('Singular matrix detected during solve.');
                    end
                catch
                    % Numerical issue: send a message and return best-effort
                    send(q, sprintf('! [Iter %d] Solver Error: Numerical instability.', iSim));
                    wlsConverged = false;
                    break;
                end
                dx = alpha * [0; dx_red];
                dVa_deg = rad2deg(dx(1:nBus));
                dVm     = dx(nBus+1:end);
                xCurrent.VA = xCurrent.VA + dVa_deg;
                xCurrent.VM = xCurrent.VM + dVm;
                max_update = max(abs(dx));
                if max_update < tol
                    wlsConverged = true;
                    convIter = iter;
                    stateOut = xCurrent;
                    break;
                end
            end
            stateOut = xCurrent;
            send(q, iSim);
        end

        function estData = estCorruptData(obj, data, Options)
            % ESTCORRUPTDATA Replaces corrupted data with a synthesized estimate.
            % Selects between 'Historical' (Phantom Anchor) and 'MME' (Hybrid Fusion)
            % based on the Options.dataImpute parameter.
          % ESTCORRUPTDATA - Replace corrupted sensor readings using GMM anchors
          %
          % Input arguments:
          % obj     - object with sensorData per bus
          % data    - table/struct with fields: Corrupt, Bus, DataType, Measurements, Variances
          % Options - struct with field dataImpute ('MME' or 'MMSE')
          %
          % Output arguments:
          % estData - same as data with corrupted entries imputed
            estData = data;
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
                    % chainName = "PF-QF-" + extractAfter(dataType, "-");
                    chainName = "PF-QF";
                    dimIdx = 1;
                elseif startsWith(dataType, "QF")
                    % chainName = "PF-QF-" + extractAfter(dataType, "-");
                    chainName = "PF-QF";
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
                gm     = chainObj.Model;
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
                % ROUTING: Select Filtering Methodology
                % ---------------------------------------------------------
                if strcmpi(Options.dataImpute, 'MME')
                    % Use historical anchor or large-noise pseudo-measurement
                    % --- ITEM 4: The "Phantom Anchor" Method ---
                    if chainName == "VM"
                        % Voltages are stable; use the tight clean variance
                        estData.Measurements(rowIdx) = mu_clean;
                        estData.Variances(rowIdx)    = var_clean;
                    else
                        % Power is dynamic; use clean mean but massive noise variance
                        % zPseudo = mu_clean / basVal;
                        zPseudo = mu_clean;
                        if dataType == "PD" || dataType == "QD"
                            zPseudo = -abs(zPseudo);
                        end
                        estData.Measurements(rowIdx) = zPseudo;
                        % Mathematically mimics Huber down-weighting
                        % varPseudo = var_noise / (basVal^2);
                        varPseudo = var_noise;
                        estData.Variances(rowIdx) = varPseudo;
                    end
                 
                elseif strcmpi(Options.dataImpute, 'MMSE')
                    % Bayesian fusion between measurement and GMM components
                    % --- ITEM 5: 1D Bayesian MMSE Fusion ---
                                      
                    if chainName == "VM"
                        % Get original measurement 
                        meas = estData.Measurements(rowIdx);
                        postProb = posterior(gm, meas);
                        postClean = postProb(cleanIdx);
                        postNoise = postProb(noiseIdx);
                        sigma2Clean = sigma2(:,:,cleanIdx);
                        sigma2noise = sigma2(:,:,noiseIdx);
                        % Voltages are stable; bypass fusion and use the clean anchor
                        estData.Measurements(rowIdx) = mu_clean;
                        % estData.Variances(rowIdx)    = var_clean;
                        estData.Variances(rowIdx)    = postClean * sigma2Clean + postNoise * sigma2noise;
                    else
                        % Build the measurement vector corresponding to P/Q pair
                        if contains(dataType, 'P')
                            meas = estData.Measurements(rowIdx:rowIdx+1)';
                        elseif contains(dataType, 'Q')
                            meas = estData.Measurements(rowIdx-1:rowIdx)';
                        end
                        postProb = posterior(gm, meas);
                        postClean = postProb(cleanIdx);
                        postNoise = postProb(noiseIdx);
                        sigma2Clean = diag(sigma2(:,:,cleanIdx));
                        sigma2noise = diag(sigma2(:,:,noiseIdx));
                        
                        covPseudo = postClean * sigma2Clean + postNoise * sigma2noise;  
                        zPseudo = mu_clean;
                        if dataType == "PD" || dataType == "QD"
                            zPseudo = -abs(zPseudo);
                        end
                        estData.Measurements(rowIdx) = zPseudo;
                        if dataType == "PD"
                            varPseudo = covPseudo(1);
                        elseif dataType == "QD"
                            varPseudo = covPseudo(2);
                        end
                     
                                           
                             
                        estData.Measurements(rowIdx) = zPseudo;
                        
                        % The "Low Confidence" Variance 
                        % Use massive noise variance to down-weight in WLS Jacobian
                       
                        % varPseudo = var_noise;
                        estData.Variances(rowIdx) = max(varPseudo, 1e-6);
                    end
                    
                else
                    error('Unknown imputation method. Set Options.dataImpute to ''MME'' or ''MMSE''.');
                end
            end
        end

        function W_iter = applyRobustWeights(~, W_base, raw_res, Options)
          % APPLYROBUSTWEIGHTS - Adjust IRWLS weights using Huber robustification
          %
          % Input arguments:
          % ~         - unused object placeholder
          % W_base    - base weight matrix (precision, typically diag(1/sigma^2))
          % raw_res   - residual vector
          % Options   - struct with fields 'RobustEstimator' and 'cHuberVal'
          %
          % Output arguments:
          % W_iter    - updated weight matrix after applying robust weights
            % APPLYROBUSTWEIGHTS Adjusts the weight matrix for IRWLS based on Zhao et al.
            % Check if robust estimation is requested
            if isfield(Options, 'RobustEstimator') && strcmpi(Options.RobustEstimator, 'Huber')
                % 1. Retrieve the Tuning Constant
                % Zhao et al. suggest 1.0 to 1.5 for impulse noise
                huber_c = Options.cHuberVal;
                % 2. Robust Scale Estimation (s)
                % The paper often uses the Median Absolute Deviation (MAD) of residuals
                % as a robust scale to standardize the residuals.
                % 1.4826 is the correction factor for Gaussian consistency.
                s_scale = 1.4826 * median(abs(raw_res)) + 1e-6;
                % 3. Standardize Residuals (r_s)
                % Note: If your W_base is diag(1/sigma^2), then sqrt(diag(W_base)) is 1/sigma.
                % We multiply by the square root of W_base to normalize the residual
                % by its expected thermal standard deviation before applying the Huber threshold.
                sigma_inv = sqrt(diag(W_base));
                r_standardized = (raw_res .* sigma_inv) / s_scale;
                % 4. Compute the Huber Weighting Matrix (Q)
                % According to the paper: q(ri) = psi(ri)/ri
                % This effectively caps the influence of any residual exceeding the threshold.
                q_weights = ones(size(r_standardized));
                outlier_idx = abs(r_standardized) > huber_c;
                % Influence bounding: weights become c / |r_std|
                q_weights(outlier_idx) = huber_c ./ abs(r_standardized(outlier_idx));
                % 5. Construct final W_iter
                % The total weight used in the normal equations is W_base * Q
                % This preserves the thermal precision while suppressing the 5% impulsive spikes.
                W_iter = W_base * diag(q_weights);
            else
                % Fallback for Ordinary WLS
                W_iter = W_base;
            end
        end
 
        % Helper functions
        function updateMCProgress(obj, data, nSim)
          % UPDATEMCPROGRESS - Update central progress counter and log messages
          %
          % Input arguments:
          % obj   - object holding simsCompleted counter
          % data  - either a log string or progress payload
          % nSim  - total number of simulations expected
            % 1. Handle string logging (Warnings/Errors) immediately
            if ischar(data)
                fprintf('%s\n', data);
                return;
            end
            % 2. Handle progress tracking
            % Increment the central counter every time ANY worker finishes
            obj.simsCompleted = obj.simsCompleted + 1;
            % Update at most 20 times (approx every 5%); ensure at least once
            updateInterval = max(1, floor(nSim / 20));
            % Check progress based on the CENTRAL counter, not the worker's ID
            if mod(obj.simsCompleted, updateInterval) == 0 || obj.simsCompleted == nSim
                fprintf('Completed %d / %d simulations (%.0f%%)\n', ...
                    obj.simsCompleted, nSim, (obj.simsCompleted/nSim)*100);
            end
        end

        function W = computeWeightMatrix(~, zTable, Options)
          % COMPUTEWEIGHTMATRIX - Build diagonal measurement weight matrix
          %
          % Input arguments:
          % ~       - unused object placeholder
          % zTable  - table of measurements with fields Variances, specVariance, DataType, Corrupt, Override
          % Options - struct with fields weightMatrix (type) and FilteringOverride (bool)
          %
          % Output arguments:
          % W       - diagonal weight matrix (rows/cols removed for excluded measurements)
            % Remove comms variance for now
            % commsVar = commVarTable.Comm_Variance;
            if strcmp(Options.weightMatrix, 'unity')
                W = diag(ones(height(zTable),1));
               
            elseif strcmp(Options.weightMatrix, 'meas-est')
                W = diag( 1./zTable.Variances );
               
            elseif strcmp(Options.weightMatrix, 'meas-specs')
                % Use theoretical per-measurement variances from specs
                % Extract the raw theoretical variances
                varArray = zTable.specVariance;
                % Create a logical mask for all measurements that are NOT Voltage Magnitude
                isNotVM = (zTable.DataType ~= "VM");
                % Enforce a minimum variance floor (1e-4) strictly on Power measurements.
                % This limits the max weight of P/Q measurements to 10,000,
                % preventing ill-conditioned matrices and weight inversion.
                varArray(isNotVM) = max(varArray(isNotVM), 1e-4);
                % Construct the diagonal weight matrix
                W = diag(1 ./ varArray);
            else
                error('Unknown weight matrix type!')
            end
            if ~(Options.FilteringOverride)
                % Exclude measurements marked corrupted XOR overridden before weighting
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
          % PREPARESENSORMEAS - Build measurement table for state estimation
          %
          % Input arguments:
          % obj        - sensor container object with Grid and sensorData
          % busMapping - table mapping sensors to buses and datatypes
          % Idx        - time/index into sensor data chains
          %
          % Output arguments:
          % zTable     - busMapping with Measurements, Variances, Corrupt, etc.
            % Inline function to split data type into its elements
            % Grid base value
            % basVal = obj.Grid.getBaseVal();
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
            % Get sensor specification standard deivation
            gridSensorNoiseStd = obj.getGridSensorNoiseStd();
            % Initialize table
            zTable = busMapping;
            % Go through all the measurements/sensors
            nMeas = height(zTable);
            for iMeas = 1:nMeas
                % Lookup bus number
                iBus = zTable.Bus(iMeas);
                % Look up measurement data type
                datatype = zTable.DataType(iMeas);
                % datatype = SensorConnector.DataType(dataTypeInd);
                % !! Refactor
                branchStr = "-"+ string(zTable.Branch(iMeas));
                branchStr(ismissing(branchStr)) = "";
                dataType = datatype + branchStr;
                % !! Refactor code
                % Find corresponding chain
                dataTypeChainIdx = contains(chainDataTypeCell{iBus}, datatype) & ...
                    contains(chainDataTypeCell{iBus}, branchStr);
                selectedChain =  obj.sensorData(iBus).Chains(dataTypeChainIdx);
                % Find normal & courrpt GMM compnents
                postProb = selectedChain.prob;
                normalGmmCmptIdx = find(postProb == max(postProb));
                corrptGmmCmptIdx = find(postProb == min(postProb));
                % Select corresponding variance index
                if table2array(selectedChain.Corrupt(Idx,:))
                    selIdx = corrptGmmCmptIdx;
                else
                    selIdx = normalGmmCmptIdx;
                end
                % Assign measurements and variances
                if strcmp(dataType, "VM")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,:);
                    zTable.Variances(iMeas) = selectedChain.Variances(:,:,selIdx);
                elseif strcmp(dataType, "PG") || strcmp(dataType, "PD") || strcmp(extractBetween(dataType,1,2), "PF")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,1);
                    zTable.Variances(iMeas) = selectedChain.Variances(1,1,selIdx);
                elseif strcmp(dataType, "QG") || strcmp(dataType, "QD") || strcmp(extractBetween(dataType,1,2), "QF")
                    zTable.Measurements(iMeas) = selectedChain.data(Idx,2);
                    zTable.Variances(iMeas) = selectedChain.Variances(2,2,selIdx);
                end
                % Assign the sensor specs variance
                zTable.specVariance(iMeas) = gridSensorNoiseStd{iBus}.(dataType).^2;
                % Assign corrupt flag
                zTable.Corrupt(iMeas) = table2array(selectedChain.Corrupt(Idx,:));
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
                    % zTable.Corrupt(pgIdx) = true; % !!
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
                    % zTable.Corrupt(qgIdx) = true; % !!
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


           pgInx = (zTable.DataType(:) == "PG");
           qgInx = (zTable.DataType(:) == "QG");

            % Get corresponding buses
            pgBus = zTable.Bus(pgInx,:);
            qgBus = zTable.Bus(qgInx,:);

            % Check if there are PD
            for iBus = 1:numel(pgBus)
                busIdx = (zTable.Bus == pgBus(iBus));

                if any(zTable.DataType(busIdx,:) == "PD")
                    % Combine if there are (PG - PD)
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
          % GETBUSTYPEIND - Return indices of buses matching a given type
            buses = {obj.Grid.Bus.Type};
            busTypeInd = find(string(buses) == busType);
        end
        function [Pinj, Qinj] = collectPQinj(obj, data)
          % COLLECTPQINJ - Extract and clean per-unit active/reactive injections
          % Input arguments:
          % data - measurement table (already per-unit)
          %
          % Output arguments:
          % Pinj - active power injections with slack removed
          % Qinj  - reactive power injections with slack and PV removed
            slackBusInd = obj.getBusTypeInd('Ref');
            pvBusInd = obj.getBusTypeInd('PV');
            Ptot = obj.extractMeas(data,"PD");
            Qtot = obj.extractMeas(data,"QD");
            % Replace unavailable values with nominal values
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
          % COLLECTPQINJOLD - Collect injected active and reactive powers from data
          %
          % Input arguments:
          % obj  - object providing grid metadata and helpers
          % data - measurement container used by extractMeas
          %
          % Output arguments:
          % Pinj - vector of active power injections with slack removed
          % Qinj  - vector of reactive power injections with slack and PV removed
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
          % PREPARECOMMVAR - Build table of communication variances per bus/type
          %
          % Input arguments:
          % obj - object containing sensorData and helper methods
          % Idx - index selecting a row from comm variance matrices
          %
          % Output arguments:
          % commVarTable - table with Bus, DataType, and Comm_Variance columns
            %
            % Inline function to split data type into its elements
            splitPair = @(s) split(s, "-");
            commVarCell = {};
            % Idx = 1; % data index
            nBus = numel(obj.sensorData);
            for iBus = 1:nBus
                nChains = numel(obj.sensorData(iBus).Chains);
                for iChains = 1:nChains
                    % Bus ID
                    busID = obj.getBusID(iBus);
                    % Chain data type
                    chainDataType = obj.getChainDataType(iBus, iChains);
                    % Comms variance
                    commVar = obj.getCommVar(iBus, iChains);
                    % if strcmp(obj.sensorData(iBus).Chains(iChains).dataType ,"VM")
                    if strcmp(chainDataType  ,"VM")
                        r = 1; % rows
                    else
                        r = 2; % rows
                    end
                    commVarCell(end + 1, :) = {busID*ones(r,1), ...
                        splitPair(chainDataType), ...
                        commVar(Idx,:)'};
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
          % COMPUTEMEAS - Predict measurements from state vector and network model
          %
          % Input arguments:
          % obj - object containing Grid with network matrices and mpc
          % x   - state struct with VM (magnitudes) and VA (angles, degrees)
          % z   - table of measurements (with Bus, Branch, DataType columns)
          %
          % Output arguments:
          % zHatTable - input table augmented with predicted Measurements and Override
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
            % DEBUG flag
            DEBUG = false;
            for iData = 1:nData
                idx = double(z{iData, 'Bus'});
                switch char(z{iData, 'DataType'})
                    case 'VM'
                        % Return voltage magnitude from state
                        zHat(iData) = x.VM(idx);
                    case 'VA'
                        % Return voltage angle from state (degrees)
                        zHat(iData) = x.VA(idx);
                    case {'PD', 'PG'}
                        % Predicted active power injection at bus
                        zHat(iData) = P_inj_all(idx);
                        % DEBUG
                        if DEBUG & (abs(zHat(iData)) > 10)
                            fprintf('%s = %1.4f\n', z{iData, 'DataType'}, zHat(iData));
                            keyboard
                        end
                    case {'QD','QG'}
                        % Predicted reactive power injection at bus
                        zHat(iData) = Q_inj_all(idx);
                        % DEBUG
                        if DEBUG & (abs(zHat(iData)) > 10)
                            fprintf('%s = %1.4f\n', z{iData, 'DataType'}, zHat(iData));
                            keyboard
                        end
                    case 'PF'
                        % Predicted active power flow on branch (from bus)
                        branchIdx = double(z{iData, 'Branch'});
                        zHat(iData) = P_flow_f_all(branchIdx);
                    case 'QF'
                        % Predicted reactive power flow on branch (from bus)
                        branchIdx = double(z{iData, 'Branch'});
                        zHat(iData) = Q_flow_f_all(branchIdx);
                end
            end
            zHatTable = z;
            zHatTable.Measurements = zHat;
            % Add Override column
            zHatTable.Override = false(nData, 1);
           
        end

        function r = computeResiduals(obj, zTable, zHatTable)
          % COMPUTERESIDUALS - Compute per-measurement residuals between tables
          %
          % Input arguments:
          % obj      - caller object (unused here)
          % zTable   - measured data table with fields Measurements, Corrupt, Override
          % zHatTable- predicted data table with at least Measurements column
          %
          % Output arguments:
          % r        - table with original first two columns and Residuals column
            nMeas = height(zTable);
            % Precallocate residual table
            r = zTable(:,1:2);
            r.Residuals = zeros(height(r),1);
            for iMeas = 1:nMeas
                % Ensure corresponding entries refer to the same measurement type
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
                % Exclude rows where Corrupt and Override differ (XOR)
                excludeFlag = xor(zTable.Corrupt, zTable.Override);
                r(excludeFlag,:) = [];
            end
        end

        function H = computeJacobian(obj, x, measMapping)
          % COMPUTEJACOBIAN - Compute measurement Jacobian for WLS state estimation
          %
          % Input arguments:
          % obj         - estimator object containing Grid and settings
          % x           - state vector struct with VM (magnitudes) and VA (degrees)
          % measMapping - table mapping measurements to buses/branches and types
          %
          % Output arguments:
          % H - sparse Jacobian matrix (measurements x [angles; magnitudes])
            % COMPUTEJACOBIAN Computes H matrix for WLS State Estimation
            
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
          % EXTRACTMEAS - Extract a specific measurement type for all buses
          %
          % Input arguments:
          % obj - object containing Grid with Bus list
          % measurements - struct/table with measurement fields
          % dataType - desired measurement type to extract
          %
          % Output arguments:
          % meas - nBus-by-1 vector of measurements (NaN where missing/corrupt)
            % Go through the measurements and extract measurement dataType
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
          % FILTEROVERRIDE - Mark measurement overrides based on corruption/settings
          %
          % If FilteringOverride is true, set override for all rows; otherwise
          % mark only fully-corrupt buses' rows as overridden.
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