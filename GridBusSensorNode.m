%%
classdef GridBusSensorNode < handle
    %GridBusSensorNode Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Bus             % struct for grid bus
        SensorConnector % struct for which data from the bus is selected
        Chains          % Processing chains
        baseVal         % Grid base value
    end

    methods
        %% Constructor
        function obj = GridBusSensorNode(Grid, iBus, dbscanGmmParam)
            % Constructor for GridBusSensorNode

            obj.baseVal = Grid.baseVal;
            obj.Bus = Grid.Bus(iBus);
            SensorConnector = Grid.sensorTable;

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
            % chainParams = struct('dbscanParams', struct('eps', 0.6, 'MinPts', 5), ...
            %                      'gmmParams',    struct('NumComponents', 2, ...
            %                                             'Replicates', 1, ...
            %                                             'RegularizationValue',1e-9, ...
            %                                             'CovarianceType', 'diagonal'));

            dataTypeCnt = 1;
            for iChain = 1:nChains

                obj.Chains(iChain).dbscanModel.Params  = chainParams.dbscanParams;
                obj.Chains(iChain).dbscanModel.Corepts = [];
                obj.Chains(iChain).dbscanModel.Idx     = [];
                obj.Chains(iChain).gmmModel.Params     = chainParams.gmmParams;

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

            iBus = obj.Bus.id;
        end

        function basVal = getBasVal(obj)

            basVal = obj.baseVal;
        end

        function meas = getMeas(obj, reqDataType, dataRange)

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
        

        % Chain parameters
        function dataType = getChainDataType(obj, ChainIdx)
            
            dataType = obj.Chains(ChainIdx).name;
        end

        % GMM parameters
        function gm = getGmmModel(obj, ChainIdx)
            
            gm = obj.Chains(ChainIdx).gmmModel.Model;
        end

        function  NumComponents = getChainGmmNumComp(obj, ChainIdx)
             
            NumComponents = obj.Chains(ChainIdx).gmmModel.Params.NumComponents;
        end

        function Replicates = getChainGmmReplicates(obj, ChainIdx)
            
            Replicates = obj.Chains(ChainIdx).gmmModel.Params.Replicates;
            
        end

        function RegularizationValue = getChainGmmRegValue(obj, ChainIdx)
            
            RegularizationValue = obj.Chains(ChainIdx).gmmModel.Params.RegularizationValue;
            
        end

        function filteredData = getChainGmmFilteredData(obj, ChainIdx)
            
            filteredData = obj.Chains(ChainIdx).gmmModel.filteredData;
        end

        function dataMean = getChainGmmDataMean(obj, ChainIdx)
            
            dataMean = obj.Chains(ChainIdx).gmmModel.Model.mu;
        end

        function dataVar = getChainGmmDataVar(obj, ChainIdx)
            
            dataVar = obj.Chains(ChainIdx).gmmModel.Model.Sigma;
        end

        function badDataIdx = getChainGmmBadDataIdx(obj, ChainIdx)
            
            badDataIdx = obj.Chains(ChainIdx).gmmModel.badDataIdx;
        end

        function ComponentProportion = getChainGmmComptProb(obj, ChainIdx)
            
            ComponentProportion = obj.Chains(ChainIdx).gmmModel.Model.ComponentProportion;
        end

        function Eps = getEpsDbscanModel(obj, ChainIdx)

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

            % Get data
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

            ChainIdx = strcmp([obj.Chains.name], chainDataType);

            % Get measurements
            meas = obj.getMeas(DataType, dataRange);

            % Normalize to per-unit
            % perUnitmeas = obj.normalizePerUnit(meas, chainDataType);

            % Override normalization
            % perUnitmeas = meas;

            % Run DBSCAN training
            % InitVals = obj.dbscanInitialization(ChainIdx, perUnitmeas);
            InitVals = obj.dbscanInitialization(ChainIdx, meas);

            % Run GMM training
            % obj.gmTrainning(ChainIdx, perUnitmeas, InitVals);
            obj.gmTrainning(ChainIdx, meas, InitVals);


        end

        function InitVals = dbscanInitialization(obj, ChainIdx, RawData)

            % Get DBSCAN parameters
            Eps = obj.getEpsDbscanModel(ChainIdx);
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
            [~, cSig] = size(Signal);

            if DEBUG
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
                fprintf('Zero covariance!\n')
                keyboard
            end

            InitVals.ComponentProportion = [sum(corePts)/length(corePts) ...
                1-sum(corePts)/length(corePts)];

            % To ensure nonzero component proportion
            InitVals.ComponentProportion(1) = InitVals.ComponentProportion(1) -eps;
            InitVals.ComponentProportion(2) = InitVals.ComponentProportion(2) +eps;          

        end

        function gmTrainning(obj, ChainIdx, RawData, InitVal)

                      
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

            nChains = numel(obj.Chains);

            for iChains = 1:nChains
                chainDataType = obj.getChainDataType(iChains);
                obj.gmChainFilter(chainDataType, dataRange);
            end
        end

        function gmChainFilter(obj, chainDataType, dataRange)

        
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

            % Normalize to per-unit
            % perUnitmeas = obj.normalizePerUnit(meas, chainDataType);
            % perUnitmeas = meas;
            
            % Filter data
            % goodDataIdx = obj.filterGMM(ChainIdx, perUnitmeas);    
            goodDataIdx = obj.filterGMM(ChainIdx, meas);
                   
            % Store bad data indices
            obj.Chains(ChainIdx).gmmModel.badDataIdx = ~goodDataIdx;
         
            % Store data
            % obj.Chains(ChainIdx).gmmModel.filteredData = perUnitmeas;
            obj.Chains(ChainIdx).gmmModel.filteredData = meas;

         
            % DEBUGGING
            % DEBUG = true;
            DEBUG = false;           

            if DEBUG
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


            % Visualize the DBSCAN model
            Signal = obj.Chains(iChains).dbscanModel.TrainningData.Signal;
            Noise = obj.Chains(iChains).dbscanModel.TrainningData.Noise;

            if strcmp(obj.Chains(iChains).name, "VM")
                scatter(figHandle, Signal, zeros(size(Signal)), 36, 'filled');
                hold on
                scatter(figHandle, Noise, zeros(size(Noise)), 36, 'filled', 'MarkerEdgeColor','w');
            else

                scatter(figHandle, Signal(:,1), Signal(:,2), 36, 'filled');
                hold on
                if ~isempty(Noise)
                    scatter(figHandle, Noise(:,1), Noise(:,2), 20, 'k', 'filled', 'MarkerEdgeColor','w');
                    hold on

                end
            end

            if ~strcmp(obj.Chains(iChains).name, "VM")
                gm = obj.Chains(iChains).gmmModel.Model;

                theta = linspace(0,2*pi,100);
                conf = chi2inv(0.999,2);                % scale for 99.9% ellipse
                % conf = chi2inv(0.95,2);                % scale for 95% ellipse
                for j = 1:gm.NumComponents
                    Sigma = gm.Sigma(:,:,j);
                    [V,D] = eig(Sigma);
                    a = sqrt(conf*diag(D));            % axis lengths
                    ellipse = (V*(a.*[cos(theta); sin(theta)])) + gm.mu(j,:)';
                    plot(figHandle, ellipse(1,:), ellipse(2,:), 'LineWidth',1.5);
                end
            end

            title([obj.Chains(iChains).name, ' Measure Data'])
            grid on
            axis(figHandle, 'equal');
        end

        function txData = transmit(obj)

            nChains = numel(obj.Chains);

            txData = struct('Bus_ID', [], 'Chains',[]);
            txData.Chains = repmat(struct('dataType', [],  ...
                                          'Means'   , [] ,  ...
                                          'Variances',[], ...
                                          'prob'     ,[], ...
                                          'data',     [], ...
                                          'Corrupt'  ,[]), 1, nChains);

            txData.Bus_ID   = obj.Bus.id;

            for iChains = 1:nChains
                % Extract data from chain getChainDataType(obj, ChainIdx)
                % txData.Chains(iChains).dataType  = obj.Chains(iChains).name;
                txData.Chains(iChains).dataType  = obj.getChainDataType(iChains);
                % txData.Chains(iChains).data      = obj.Chains(iChains).gmmModel.filteredData;
                txData.Chains(iChains).data      = obj.getChainGmmFilteredData(iChains);
                % txData.Chains(iChains).Means     = obj.Chains(iChains).gmmModel.Model.mu;
                txData.Chains(iChains).Means     = obj.getChainGmmDataMean(iChains);
                % txData.Chains(iChains).Variances = obj.Chains(iChains).gmmModel.Model.Sigma;
                txData.Chains(iChains).Variances = obj.getChainGmmDataVar(iChains);
                % txData.Chains(iChains).prob      = obj.Chains(iChains).gmmModel.Model.ComponentProportion;
                txData.Chains(iChains).prob      = obj.getChainGmmComptProb(iChains);
                % txData.Chains(iChains).Corrupt   = obj.Chains(iChains).gmmModel.badDataIdx;
                txData.Chains(iChains).Corrupt   = obj.getChainGmmBadDataIdx(iChains);

            end
        end

        %% Helper methods
        function perUnitmeas = normalizePerUnit(obj, meas, chainDataType)
          
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