%%
classdef GridBusSensorNode < handle
    %GRIDBUSSENSORNODE - Representation of a sensor node attached to a grid bus
    %   OBJ = GRIDBUSSENSORNODE() encapsulates the bus data, selection
    %   connector, processing chains, and base values used by a sensor node.
    %
    %   Properties:
    %       Bus             - struct containing grid bus data
    %       SensorConnector - struct defining which bus data to select
    %       Chains          - processing chains applied to sensor data
    %       baseVal         - grid base value for normalization/scaling
    %       sensorImpNoisInd- index or flag for imperfect/noisy sensor handling
    %   See also HANDLE
    %GridBusSensorNode Summary of this class goes here
    %   Detailed explanation goes here
    properties
        Bus             % struct for grid bus
        SensorConnector % struct for which data from the bus is selected
        Chains          % Processing chains
        baseVal         % Grid base value
        sensorImpNoisInd
    end

    methods
        %% Constructor
        function obj = GridBusSensorNode(Grid, iBus, dbscanGmmParam)
            % GRIDBUSSENSORNODE - Construct sensor node for a specific bus
            %
            % Input arguments:
            % Grid - grid struct containing bus and sensor metadata
            % iBus - index of the bus to attach sensors to
            % dbscanGmmParam - struct of parameters for dbscan and GMM models
            % Constructor for GridBusSensorNode
            obj.baseVal = Grid.baseVal;
            obj.Bus = Grid.Bus(iBus);
            SensorConnector = Grid.sensorTable;
            obj.sensorImpNoisInd = [];
            % Convert table to struct
            busNum = obj.Bus.id; % Bus number
            % Get available sensor data types
            dataTypeInd = ( SensorConnector.Bus == busNum);
            % dataType = SensorConnector.DataType(dataTypeInd);
            % !! Refactor
            datatype = SensorConnector.DataType(dataTypeInd);
            branchStr = "_"+ string(SensorConnector.Branch(dataTypeInd));
            branchStr(ismissing(branchStr)) = "";
            dataType = datatype + branchStr;
            % Check if VM is there
            vmInd = double(any(ismember(dataType, "VM")));
            % Compute the number of chains
            nChains = (numel(dataType)+ vmInd)/2 ;
            % Chain parameters
            chainParams = dbscanGmmParam;
            dataTypeCnt = 1;
            for iChain = 1:nChains
                obj.Chains(iChain).dbscanModel.Params  = chainParams.dbscanParams;
                obj.Chains(iChain).dbscanModel.Corepts = [];
                obj.Chains(iChain).dbscanModel.Idx     = [];
                obj.Chains(iChain).gmmModel.Params     = chainParams.gmmParams;
                obj.Chains(iChain).badDataIdx          = [];
                % Assign chain type based on next available dataType token
                if strcmp(dataType(dataTypeCnt), "VM")
                    obj.Chains(iChain).name = "VM";
                    dataTypeCnt = dataTypeCnt + 1;
                elseif  strcmp(dataType(dataTypeCnt), "PD") || strcmp(dataType(dataTypeCnt), "QD")
                    obj.Chains(iChain).name = "PD-QD";
                    dataTypeCnt = dataTypeCnt + 2;
                elseif  strcmp(dataType(dataTypeCnt), "PG") || strcmp(dataType(dataTypeCnt), "QG")
                    obj.Chains(iChain).name = "PG-QG";
                    dataTypeCnt = dataTypeCnt + 2;
                    % Might have several PF-QF chains
                elseif  strcmp(extractBefore(dataType(dataTypeCnt),3), "PF") || strcmp(extractBefore(dataType(dataTypeCnt),3), "QF")
                    obj.Chains(iChain).name = "PF-QF-"+ extractAfter(dataType(dataTypeCnt),"_");
                    dataTypeCnt = dataTypeCnt + 2;
                end
            end
        end

        %% Getters
        function  iBus = getBusID(obj)
            % GETBUSID - Return the numeric bus identifier
            iBus = obj.Bus.id;
        end
        function basVal = getBasVal(obj)
            % GETBASVAL - Return the stored base value for this object
            basVal = obj.baseVal;
        end
        function meas = getMeas(obj, reqDataType, dataRange)
            % GETMEAS - Retrieve requested measurement columns for given rows
            %   reqDataType - string or string array of column names to fetch
            %   dataRange   - row indices or logical mask selecting rows
            % Get true values table
            BusMeas = obj.Bus.meas;
            % Get available measurements data types
            avalDataTypes = string(BusMeas.Properties.VariableNames);
            if all(ismember(reqDataType, avalDataTypes))
                meas =  BusMeas{dataRange,reqDataType};
            else
                unAvalDataTypeIdx = ~ismember(reqDataType, avalDataTypes);
                unAvalDataType = reqDataType(unAvalDataTypeIdx);
                fprintf('Measurement data type %s is/are not available at bus %d\n', char(strjoin(unAvalDataType, ', ')), iBus);
                fprintf('Available data is\n');
                disp(BusMeas);
                error('Measurement data type mismatch!')
            end
        end

        function sensorImpNoisInd = getSensorImpNoisInd(obj)
            % GETSENSORIMPNOISIND - Return sensor impedance/noise indicator property
            %
            % Output:
            % sensorImpNoisInd - stored sensor impedance/noise indicator
            sensorImpNoisInd = obj.sensorImpNoisInd;
        end
        % Chain parameters
        function dataType = getChainDataType(obj, ChainIdx)
            % GETCHAINDATATYPE - Return the name/type for the specified chain
            %
            % Inputs:
            % ChainIdx - index of chain to query
            %
            % Output:
            % dataType - chain name/type
            dataType = obj.Chains(ChainIdx).name;
        end
        function badDataIdx = getBadDataIdx(obj, ChainIdx)
            % GETBADDATAIDX - Return indices of bad data for specified chain
            %
            % Inputs:
            % ChainIdx - index of chain to query
            %
            % Output:
            % badDataIdx - indices flagged as bad
            badDataIdx = obj.Chains(ChainIdx).badDataIdx;
        end

        % GMM parameters
        function gm = getGmmModel(obj, ChainIdx)
            % GETGMMMODEL - Return the fitted GMM object for a chain
            %
            % Input arguments:
            % obj      - container with Chains
            % ChainIdx - index of the chain
            gm = obj.Chains(ChainIdx).gmmModel.Model;
        end
        function  NumComponents = getChainGmmNumComp(obj, ChainIdx)
            % GETCHAINGMMNUMCOMP - Number of mixture components in chain GMM
            NumComponents = obj.Chains(ChainIdx).gmmModel.Params.NumComponents;
        end
        function Replicates = getChainGmmReplicates(obj, ChainIdx)
            % GETCHAINGMMREPLICATES - Number of replicates used to fit chain GMM
            Replicates = obj.Chains(ChainIdx).gmmModel.Params.Replicates;

        end

        function RegularizationValue = getChainGmmRegValue(obj, ChainIdx)
            % GETCHAINGMMREGVALUE - Return GMM regularization value for a chain
            %
            % Input arguments:
            % ChainIdx - index of the chain to query
            %
            % Output arguments:
            % RegularizationValue - stored regularization hyperparameter
            RegularizationValue = obj.Chains(ChainIdx).gmmModel.Params.RegularizationValue;

        end
        function Model = getChainModel(obj, ChainIdx)
            % GETCHAINMODEL - Return the fitted GMM model for a chain
            %
            % Model - structure/object representing the GMM fit
            Model = obj.Chains(ChainIdx).gmmModel.Model;
        end
        function filteredData = getChainGmmFilteredData(obj, ChainIdx)
            % GETCHAINGMMFILTEREDDATA - Retrieve data after GMM preprocessing/filtering
            filteredData = obj.Chains(ChainIdx).gmmModel.filteredData;
        end
        function dataMean = getChainGmmDataMean(obj, ChainIdx)
            % GETCHAINGMMDATAMEAN - Return component means (mu) of the chain's GMM
            dataMean = obj.Chains(ChainIdx).gmmModel.Model.mu;
        end
        function dataVar = getChainGmmDataVar(obj, ChainIdx)
            % GETCHAINGMMDATAVAR - Return component covariances (Sigma) of the chain's GMM
            dataVar = obj.Chains(ChainIdx).gmmModel.Model.Sigma;
        end

        function badDataIdx = getChainGmmBadDataIdx(obj, ChainIdx)
            % GETCHAINGMMBADDATAIDX - Return indices of bad data for a chain's GMM
            %
            % Input arguments:
            % obj      - object containing Chains
            % ChainIdx - index of the chain
            %
            % Output arguments:
            % badDataIdx - indices of data marked as bad for that chain
            badDataIdx = obj.Chains(ChainIdx).badDataIdx;
        end
        function ComponentProportion = getChainGmmComptProb(obj, ChainIdx)
            % GETCHAINGMMCOMPTPROB - Get GMM component proportions for a chain
            %
            % Input arguments:
            % obj      - object containing Chains
            % ChainIdx - index of the chain
            %
            % Output arguments:
            % ComponentProportion - vector of GMM component weights
            ComponentProportion = obj.Chains(ChainIdx).gmmModel.Model.ComponentProportion;
        end
        function Eps = getEpsDbscanModel(obj, ChainIdx)
            % GETEPSDBSCANMODEL - Return epsilon parameter from a chain's DBSCAN
            %
            % Input arguments:
            % obj      - object containing Chains
            % ChainIdx - index of the chain
            %
            % Output arguments:
            % Eps - DBSCAN epsilon neighborhood radius
            Eps = obj.Chains(ChainIdx).dbscanModel.Params.eps;
        end

        function MinPts = getMinPtsDbscanModel(obj, ChainIdx)

            MinPts = obj.Chains(ChainIdx).dbscanModel.Params.MinPts;
        end
        %%
        function trainChains(obj, dataRange)
            % TRAINCHAINS - Train all chains with provided training count
            %
            % Input arguments:
            % obj - object containing Chains array
            % nTrainningData - number of training samples per chain

            for iChain = 1:numel(obj.Chains)
                chainDataType = obj.getChainDataType(iChain);
                obj.trainChain(chainDataType, dataRange);
            end
        end

        function trainChain(obj, chainDataType, dataRange)
            % TRAINCHAIN - Train chain models for a specified data type and range
            %
            % Input arguments:
            % chainDataType - identifier for chain measurement types (string)
            % dataRange     - time or index range for measurements
            % Get data
            if strcmp(chainDataType, "Ref") || strcmp(chainDataType, "VM")
                DataType = {'VM'};
            elseif strcmp(chainDataType,"PD-QD")
                DataType = {'PD', 'QD'};
            elseif strcmp(chainDataType,"PG-QG")
                DataType = {'PG', 'QG'};
            elseif contains(chainDataType,"PF-QF")
                % Extract brand ID after second hyphen to build PF/QF types
                hyphenInd = strfind(chainDataType, "-");
                Name = char(chainDataType);
                brandID = Name(hyphenInd(2)+1:end);
                DataType = {char('PF-'+string(brandID)), char('QF-'+string(brandID))};
            end
            ChainIdx = strcmp([obj.Chains.name], chainDataType);
            % Get measurements
            meas = obj.getMeas(DataType, dataRange);

            % Run DBSCAN training
            InitVals = obj.dbscanInitialization(ChainIdx, meas);
            % Run GMM training
            obj.gmTrainning(ChainIdx, meas, InitVals);
        end

        function InitVals = dbscanInitialization(obj, ChainIdx, RawData)
            % DBSCANINITIALIZATION - Initialize GMM parameters using DBSCAN clustering
            %
            % Input arguments:
            % obj      - object containing chain configurations
            % ChainIdx - index of the chain to initialize
            % RawData  - observations to cluster (rows = samples)
            %
            % Output arguments:
            % InitVals - struct with fields mu, Sigma, ComponentProportion
            % Get DBSCAN parameters
            Eps    = obj.getEpsDbscanModel(ChainIdx);
            MinPts = obj.getMinPtsDbscanModel(ChainIdx);

            % Z-scoring
            normData = zscore(RawData);
            % DBSCAN filtering
            [Idx, corePts] = dbscan(normData, Eps, MinPts);
            % Store DBSCAN output
            obj.Chains(ChainIdx).dbscanModel.Corepts = corePts;
            obj.Chains(ChainIdx).dbscanModel.Idx = Idx;
            % Filtered data
            Signal = RawData(Idx == 1,:);
            Noise = RawData(Idx ~= 1,:);
            % Store trainning data
            obj.Chains(ChainIdx).dbscanModel.TrainningData.Signal = Signal;
            obj.Chains(ChainIdx).dbscanModel.TrainningData.Noise = Noise;
            % DEBUGING
            % DEBUG = true;
            DEBUG = false;
            if DEBUG
                [~, cSig] = size(Signal);
                figure
                if cSig == 1
                    scatter(Signal, zeros(size(Signal)))
                    hold on, grid on
                    scatter(Noise, zeros(size(Noise)))
                else
                    scatter(Signal(:,1), Signal(:,2))
                    hold on, grid on
                    scatter(Noise(:,1), Noise(:,2))
                end
            end
            % Compute initial values
            InitVals.mu(1,:) = mean(Signal);
            InitVals.mu(2,:) = mean(Noise);
            InitVals.Sigma(:,:,1) = cov(Signal);
            InitVals.Sigma(:,:,2) = cov(Noise);
            % DEBUG
            if any(cov(Noise) == 0)
                fprintf('Zero covariance in GMM fitting!\n')
                keyboard
            end
            InitVals.ComponentProportion = [sum(corePts)/length(corePts) ...
                1-sum(corePts)/length(corePts)];
            % To ensure nonzero component proportion
            InitVals.ComponentProportion(1) = InitVals.ComponentProportion(1) -eps;
            InitVals.ComponentProportion(2) = InitVals.ComponentProportion(2) +eps;
        end

        function gmTrainning(obj, ChainIdx, RawData, InitVal)
            % GMTRAINNING - Fit and store a Gaussian Mixture Model for a chain
            %
            % Input arguments:
            % obj      - object containing chain configurations and storage
            % ChainIdx - index of the chain to train
            % RawData  - N-by-D data matrix used to fit the GMM
            % InitVal  - initialization structure/values for fitgmdist

            % GMM parameters
            NumComponents = obj.getChainGmmNumComp(ChainIdx);

            Replicates = obj.getChainGmmReplicates(ChainIdx);
            RegularizationValue = obj.getChainGmmRegValue(ChainIdx);
            % GMM training
            opts = statset('MaxIter',1000,'TolFun',1e-5);
            gm = fitgmdist(RawData, NumComponents, ...
                'Start', InitVal, ...
                'Replicates', Replicates, ...
                'RegularizationValue',RegularizationValue, ...
                'CovarianceType', 'full', ...
                'Options', opts);
            % Store GMM model
            obj.Chains(ChainIdx).gmmModel.Model = gm;
            obj.Chains(ChainIdx).gmmModel.Data = RawData;

        end

        function gmFiltering(obj, dataRange)
            % GMFILTERING - Filter each chain's data with Gaussian mixture model
            %
            % Input arguments:
            % dataRange - time or index range used for chain filtering
            nChains = numel(obj.Chains);
            sensorImpInd = table();
            for iChains = 1:nChains
                chainDataType = obj.getChainDataType(iChains);
                obj.gmChainFilter(chainDataType, dataRange);
                chainImpInd = obj.getChainGmmBadDataIdx(iChains);
                % Combine per-chain bad-data indicator tables into one table
                if isempty(sensorImpInd.Properties.VariableNames)
                    sensorImpInd = chainImpInd;
                else
                    sensorImpInd = [sensorImpInd, chainImpInd];
                end


            end
            obj.sensorImpNoisInd = sensorImpInd;
        end

        function gmChainFilter(obj, chainDataType, dataRange)
            % GMCHAINFILTER - Filter measurements for a chain using GMM-based mask
            %
            % Input arguments:
            % obj           - object containing Chains and methods
            % chainDataType - chain identifier or compound type string
            % dataRange     - time/index range for measurement retrieval

            % Get data type
            if strcmp(chainDataType, "Ref") || strcmp(chainDataType, "VM")
                DataType = {'VM'};
            elseif strcmp(chainDataType,"PD-QD")
                DataType = {'PD', 'QD'};
            elseif strcmp(chainDataType,"PG-QG")
                DataType = {'PG', 'QG'};
            elseif contains(chainDataType,"PF-QF")
                hyphenInd = strfind(chainDataType, "-");
                Name = char(chainDataType);
                brandID = Name(hyphenInd(2)+1:end);
                DataType = {char('PF-'+string(brandID)), char('QF-'+string(brandID))};
            end
            % Get chain index
            ChainIdx = strcmp([obj.Chains.name], chainDataType);
            % Get measurements
            meas = obj.getMeas(DataType, dataRange);

            % Filter data
            goodDataIdx = obj.filterGMM(ChainIdx, meas);

            % Store bad data indices
            obj.Chains(ChainIdx).gmmModel.badDataIdx = ~goodDataIdx;
            obj.Chains(ChainIdx).badDataIdx = array2table(~goodDataIdx, ...
                'VariableNames', {char(chainDataType)});

            % Store data
            obj.Chains(ChainIdx).gmmModel.filteredData = meas;

            % DEBUGGING
            % DEBUG = true;
            DEBUG = false;
            if DEBUG
                % Split measurements for visualization into signal and noise sets
                Signal = perUnitmeas(goodDataIdx,:);
                Noise = perUnitmeas(~goodDataIdx,:);
                [~, cSig] = size(Signal);
                figure
                if cSig == 1
                    scatter(Signal, zeros(size(Signal)))
                    hold on, grid on
                    scatter(Noise, zeros(size(Noise)))
                else
                    scatter(Signal(:,1), Signal(:,2))
                    hold on, grid on
                    scatter(Noise(:,1), Noise(:,2))
                end
            end
        end

        function   goodDataIdx = filterGMM(obj, ChainIdx, meas)
            % FILTERGMM - Identify measurements assigned to the dominant GMM component
            %
            % Input arguments:
            % obj     - object with GMM models
            % ChainIdx- index selecting which GMM to use
            % meas    - N-by-D measurement matrix
            %
            % Output arguments:
            % goodDataIdx - N-by-1 logical, true when measurement's component
            %               matches model's dominant component

            % Get GMM
            gm = obj.getGmmModel(ChainIdx);
            % Find which component is good and corrupt
            [~, goodCompInd] = max(gm.ComponentProportion,[],2);
            % Compute posterior of data
            post = posterior(gm, meas);
            % Find to which component the data belongs to
            [~, goodDataColIdx] = max(post,[],2);
            goodDataIdx = (goodCompInd == goodDataColIdx);
        end

        function visualizeModels(obj)
            % VISUALIZEMODELS - Plot each chain's model in a tiled figure
            nChains = numel(obj.Chains);
            figure;
            % Preallocate 2x2 figure handles
            figHandle = gobjects(4,1);

            for iChains = 1:nChains
                figHandle(iChains) = nexttile;

                obj.visualizeChainModel(iChains, figHandle(iChains));
            end

            % Figure title
            iBus = obj.getBusID();
            sgtitle(['Bus-' num2str(iBus)])
        end


        function visualizeChainModel(obj, iChains, figHandle)
            % VISUALIZECHAINMODEL - Plot DBSCAN clusters and GMM confidence ellipses
            %
            % Input arguments:
            % obj       - object containing Chains with models
            % iChains   - index of the chain to visualize
            % figHandle - axes or figure handle to plot into
            % Visualize the DBSCAN model
            Signal = obj.Chains(iChains).dbscanModel.TrainningData.Signal;
            Noise = obj.Chains(iChains).dbscanModel.TrainningData.Noise;
            if strcmp(obj.Chains(iChains).name, "VM")
                hold on
                hScatterCore = scatter(figHandle, Signal, zeros(size(Signal)), 36, 'filled');
                % scatter(figHandle, Noise, zeros(size(Noise)), 36, 'filled', 'MarkerEdgeColor','w');

                scatter(figHandle, Noise, zeros(size(Noise)), 36, 'x', 'MarkerEdgeColor', 'w', 'LineWidth', 1.2);
                legend('Core', 'Noise', 'FontSize',12)
            else
                % 2D case: plot core points and noise separately
                scatter(figHandle, Signal(:,1), Signal(:,2), 36, 'filled');
                hold on
                if ~isempty(Noise)
                    % scatter(figHandle, Noise(:,1), Noise(:,2), 20, 'k', 'filled', 'MarkerEdgeColor','w');
                    hScatterNoise = scatter(figHandle, Noise(:,1), Noise(:,2), 20, 'x', 'MarkerEdgeColor','w', 'LineWidth', 1.2);
                    hold on
                end

            end
            if ~strcmp(obj.Chains(iChains).name, "VM")
                % Overlay GMM confidence ellipses for each mixture component
                gm = obj.Chains(iChains).gmmModel.Model;
                theta = linspace(0,2*pi,100);
                conf = chi2inv(0.990,2);
                % conf = chi2inv(0.999,2);                % scale for 99.9% ellipse
                % conf = chi2inv(0.95,2);                % scale for 95% ellipse
                for j = 1:gm.NumComponents
                    Sigma = gm.Sigma(:,:,j);
                    [V,D] = eig(Sigma);
                    a = sqrt(conf*diag(D));            % axis lengths
                    ellipse = (V*(a.*[cos(theta); sin(theta)])) + gm.mu(j,:)';
                    hConf{j} = plot(figHandle, ellipse(1,:), ellipse(2,:), 'r', 'LineWidth',1.5);
                    legend(hConf{j}, 'off')
                end

            end
            title([obj.Chains(iChains).name, ' Measure Data'])
            grid on
            axis(figHandle, 'equal');

            legend('Core', 'Noise', 'FontSize',12)
        end

        function txData = transmit(obj)
            % TRANSMIT - Package per-chain GMM info for transmission
            %
            % Input arguments:
            % obj - object containing Chains and Bus information
            %
            % Output arguments:
            % txData - struct with Bus_ID and per-chain fields (Model, Means, etc.)
            nChains = numel(obj.Chains);
            txData = struct('Bus_ID', [], 'Chains',[]);
            txData.Chains = repmat(struct('dataType', [], ...
                'Model',    [], ...
                'Means'   , [], ...
                'Variances',[], ...
                'prob'     ,[], ...
                'data',     [], ...
                'Corrupt'  ,[]), 1, nChains);
            txData.Bus_ID   = obj.Bus.id;
            for iChains = 1:nChains
                % Populate chain-level metadata and GMM-derived fields
                % Extract data from chain getChainDataType(obj, ChainIdx)
                txData.Chains(iChains).dataType  = obj.getChainDataType(iChains);
                txData.Chains(iChains).Model     = obj.getChainModel(iChains);
                txData.Chains(iChains).data      = obj.getChainGmmFilteredData(iChains);
                txData.Chains(iChains).Means     = obj.getChainGmmDataMean(iChains);
                txData.Chains(iChains).Variances = obj.getChainGmmDataVar(iChains);
                txData.Chains(iChains).prob      = obj.getChainGmmComptProb(iChains);
                txData.Chains(iChains).Corrupt   = obj.getChainGmmBadDataIdx(iChains);
            end
        end

        %% Helper methods
        function perUnitmeas = normalizePerUnit(obj, meas, chainDataType)
            % NORMALIZEPERUNIT - Convert measurements to per-unit using object's base value
            %
            % Input arguments:
            % obj           - object providing getBasVal method
            % meas          - measurement (scalar or array) to normalize
            % chainDataType - string tag; "VM" bypasses normalization
            %
            % Output arguments:
            % perUnitmeas   - measurement in per-unit (or unchanged for "VM")

            % Grid base value
            basVal = obj.getBasVal();
            % Normalize measurements to per-unit
            if strcmp(chainDataType, "VM")
                perUnitmeas = meas;
            else
                perUnitmeas = meas / basVal;
            end
        end

    end
end