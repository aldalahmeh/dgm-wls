
%%
classdef GridBusSensorNode < handle
    %GridBusSensorNode Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Bus             % struct for grid bus
        SensorConnector % struct for which data from the bus is selected
        Chains          % Processing chains
    end

    methods
        function obj = GridBusSensorNode(Grid, iBus)
            % Constructor for GridBusSensorNode

            obj.Bus = Grid.Bus(iBus);
            SensorConnector = Grid.sensorTable;
            
            % Convert table to struct
            busNum = obj.Bus.id; % Bus number
            % Get available sensor data types
            dataTypeInd = ( SensorConnector.Bus == busNum);
            % dataType = SensorConnector.DataType(dataTypeInd); 

            datatype = SensorConnector.DataType(dataTypeInd);
            branchStr = "_"+ string(SensorConnector.Branch(dataTypeInd)); 
            branchStr(ismissing(branchStr)) = "";

            dataType = datatype + branchStr; 
            
            % Check if VM is there
            vmInd = double(any(ismember(dataType, "VM")));
            % Compute the number of chains
            nChains = (numel(dataType)+ vmInd)/2 ;

            % Chain parameters
            chainParams = struct('dbscanParams', struct('eps', 0.6, 'MinPts', 5), ...
                'gmmParams',    struct('NumComponents', 2, ...
                'Replicates', 1, ...
                'RegularizationValue',1e-6));

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

            % for k = 1:numel(dataType)
            %     SensorConnectorStruct.(dataType(k)) = false;     % assign empty array to each field
            % end

            % Go through the availabe sensor data types in dataType in set
            % the corressponding field to true
            % for i = 1:length(dataType)
            %     SensorConnectorStruct.(dataType(i)) = true;
            % end
            % 
            % % obj.SensorConnector = SensorConnector;
            % obj.SensorConnector = SensorConnectorStruct;
            % 
            % chainParams = struct('dbscanParams', struct('eps', 0.6, 'MinPts', 5), ...
            %                      'gmmParams',    struct('NumComponents', 2, ...
            %                                             'Replicates', 1, ...
            %                                             'RegularizationValue',1e-6));
            % 
            % % Get Data Selection field names
            % fn = fieldnames(obj.SensorConnector);             % cell array of field names
            % vals = cellfun(@(f) obj.SensorConnector.(f), fn); % logical array (same order as fn)
            % ConnectedSensors = fn(vals);                      % cell array of names that are true
            % 
            % % Populate processing chains
            % obj.Chains = struct([]);
            % 
            % ChainCnt = 1;
            % if isempty(ConnectedSensors)
            %     obj.Chains = struct([]);   % empty struct array if none enabled
            % else
            %     if any(strcmp(ConnectedSensors, "VM"))
            %         obj.Chains(ChainCnt).name = "VM";
            %         obj.Chains(ChainCnt).dbscanModel.Params  = chainParams.dbscanParams;
            %         obj.Chains(ChainCnt).dbscanModel.Corepts = [];
            %         obj.Chains(ChainCnt).dbscanModel.Idx     = [];
            %         obj.Chains(ChainCnt).gmmModel.Params     = chainParams.gmmParams;
            % 
            %         % Number of GMM components must be 1
            %         % obj.Chains(ChainCnt).gmmModel.Params.NumComponents = 1;
            % 
            %         ChainCnt = ChainCnt + 1;
            %     end
            % 
            %     if any(strcmp(ConnectedSensors, "PD")) && any(strcmp(ConnectedSensors, "QD"))
            %         obj.Chains(ChainCnt).name = "PD-QD";
            %         obj.Chains(ChainCnt).dbscanModel.Params  = chainParams.dbscanParams;
            %         obj.Chains(ChainCnt).dbscanModel.Corepts = [];
            %         obj.Chains(ChainCnt).dbscanModel.Idx     = [];
            %         obj.Chains(ChainCnt).gmmModel.Params     = chainParams.gmmParams;
            % 
            %         ChainCnt = ChainCnt + 1;
            %     end
            % 
            %     if any(strcmp(ConnectedSensors, "PF")) && any(strcmp(ConnectedSensors, "QF"))
            %         obj.Chains(ChainCnt).name = "PF-QF";
            %         obj.Chains(ChainCnt).dbscanModel.Params  = chainParams.dbscanParams;
            %         obj.Chains(ChainCnt).dbscanModel.Corepts = [];
            %         obj.Chains(ChainCnt).dbscanModel.Idx     = [];
            %         obj.Chains(ChainCnt).gmmModel.Params     = chainParams.gmmParams;
            % 
            %         ChainCnt = ChainCnt + 1;
            %     end
            % 
            %     if any(strcmp(ConnectedSensors, "PG")) && any(strcmp(ConnectedSensors, "QG"))
            %         obj.Chains(ChainCnt).name = "PG-QG";
            %         obj.Chains(ChainCnt).dbscanModel.Params  = chainParams.dbscanParams;
            %         obj.Chains(ChainCnt).dbscanModel.Corepts = [];
            %         obj.Chains(ChainCnt).dbscanModel.Idx     = [];
            %         obj.Chains(ChainCnt).gmmModel.Data       = []; 
            %         obj.Chains(ChainCnt).gmmModel.Model      = [];
            %         obj.Chains(ChainCnt).gmmModel.Params     = chainParams.gmmParams;
            % 
            %         ChainCnt = ChainCnt + 1;
            %     end
            % 
            % end

        end

        function trainChains(obj, dataRange)
          % TRAINCHAINS - Train all chains with provided training count
          %
          % Input arguments:
          % obj - object containing Chains array
          % nTrainningData - number of training samples per chain
            for i = 1:numel(obj.Chains)
                obj.trainChain(obj.Chains(i).name, dataRange);
            end
        end

        function trainChain(obj, name, dataRange)

            % Get data
            % if strcmp(name, "V")
            if strcmp(name, "Ref") || strcmp(name, "VM")
                DataType = {'VM'};
            elseif strcmp(name,"PD-QD")
                DataType = {'PD', 'QD'};            
            elseif strcmp(name,"PG-QG")
                DataType = {'PG', 'QG'};
            elseif contains(name,"PF-QF")
                hyphenInd = strfind(name, "-");
                Name = char(name);
                brandID = Name(hyphenInd(2)+1:end);
                % !!
                DataType = {char('PF-'+string(brandID)), char('QF-'+string(brandID))};
            end
           
            ChainIdx = strcmp([obj.Chains.name], name);            
            data = obj.Bus.meas{dataRange, DataType};

            % Run DBSCAN training
            InitVals = obj.dbscanInitialization(ChainIdx, data);

            % Run GMM training
            obj.gmTrainning(ChainIdx, data, InitVals);


        end

        function InitVals = dbscanInitialization(obj, ChainIdx, RawData)

            % Z-scoring
            normData = zscore(RawData);

            % DBSCAN filtering
            [Idx, corePts] = dbscan(normData, obj.Chains(ChainIdx).dbscanModel.Params.eps, ...
                                              obj.Chains(ChainIdx).dbscanModel.Params.MinPts);
            obj.Chains(ChainIdx).dbscanModel.Corepts = corePts;
            obj.Chains(ChainIdx).dbscanModel.Idx = Idx;

            % Filtered data
            Signal = RawData(Idx == 1,:);
            Noise = RawData(Idx ~= 1,:);
            % Signal = RawData(corePts,:);
            % Noise = RawData(~corePts,:);

            % Store trainning data
            obj.Chains(ChainIdx).dbscanModel.TrainningData.Signal = Signal;
            obj.Chains(ChainIdx).dbscanModel.TrainningData.Noise = Noise;

            % Compute initial values
            InitVals.mu(1,:) = mean(Signal);
            InitVals.mu(2,:) = mean(Noise);

            InitVals.Sigma(:,:,1) = diag(var(Signal));
            InitVals.Sigma(:,:,2) = diag(var(Noise));

            InitVals.ComponentProportion = [sum(corePts)/length(corePts) ...
                1-sum(corePts)/length(corePts)];

            % To ensure nonzero component proportion
            InitVals.ComponentProportion(1) = InitVals.ComponentProportion(1) -eps;
            InitVals.ComponentProportion(2) = InitVals.ComponentProportion(2) +eps;

            % For VM data only retain one set of values
            % if strcmp(obj.Chains(ChainIdx).name, "VM")
            %     InitVals.mu = InitVals.mu(1);
            %     InitVals.Sigma = InitVals.Sigma(:,:,1);
            %     InitVals.ComponentProportion = 1;
            % end

        end

        function gmTrainning(obj, ChainIdx, RawData, InitVal)

            NumComponents = obj.Chains(ChainIdx).gmmModel.Params.NumComponents;
            Replicates = obj.Chains(ChainIdx).gmmModel.Params.Replicates;
            RegularizationValue = obj.Chains(ChainIdx).gmmModel.Params.RegularizationValue;

            gm = fitgmdist(RawData, NumComponents, 'Start', InitVal, ...
                'Replicates', Replicates, 'RegularizationValue',RegularizationValue);

            obj.Chains(ChainIdx).gmmModel.Model = gm;
            obj.Chains(ChainIdx).gmmModel.Data = RawData;

        end

        function gmFiltering(obj, dataRange)

            nChains = numel(obj.Chains);

            for iChains = 1:nChains
                obj.gmChainFilter(obj.Chains(iChains).name, dataRange);
            end
        end

        function gmChainFilter(obj, name, dataRange)
            
            % Get data           
            if strcmp(name, "Ref") || strcmp(name, "VM")
                DataType = {'VM'};
            elseif strcmp(name,"PD-QD")
                DataType = {'PD', 'QD'};
            elseif strcmp(name,"PG-QG")
                DataType = {'PG', 'QG'};
            elseif contains(name,"PF-QF")
                hyphenInd = strfind(name, "-");
                Name = char(name);
                brandID = Name(hyphenInd(2)+1:end);
                DataType = {char('PF-'+string(brandID)), char('QF-'+string(brandID))};
                % DataType = {'PF', 'QF'};
            end

            ChainIdx = strcmp([obj.Chains.name], name);

            data = obj.Bus.meas{dataRange, DataType};

            % Filter data
            gm = obj.Chains(ChainIdx).gmmModel.Model;
            post = posterior(gm, data);
            goodDataIdx = logical(post(:,1));
            % filteredData = data( goodDataIdx ,:);
            % Send NaNs for corrupted data
            % badDataVec = double(~goodDataIdx);
            % badDataVec(~goodDataIdx) = nan;

            % Store bad data indices
            obj.Chains(ChainIdx).gmmModel.badDataIdx = ~goodDataIdx;
            % Store filtered data
            obj.Chains(ChainIdx).gmmModel.filteredData = data;
            % obj.Chains(ChainIdx).gmmModel.filteredData = filteredData;
            % obj.Chains(ChainIdx).gmmModel.filteredData  ...
            %     = data + badDataVec;

        end

        function visualizeModels(obj)

            nChains = numel(obj.Chains);

            % Preallocate 2x2 figure handles
            f = figure;
            % t = tiledlayout(f, 2,2);
            figHandle = gobjects(4,1);

            sgtitle(['Bus-', string(obj.Bus.id)])

            for iChains = 1:nChains
                figHandle(iChains) = nexttile;
                obj.visualizeChainModel(iChains, figHandle(iChains));
            end
        end

        function visualizeChainModel(obj, iChains, figHandle)


            % Visualize the DBSCAN model
            Signal = obj.Chains(iChains).dbscanModel.TrainningData.Signal;
            Noise = obj.Chains(iChains).dbscanModel.TrainningData.Noise;

            if strcmp(obj.Chains(iChains).name, "V")
                scatter(figHandle, Signal, zeros(size(Signal)), 36, 'filled');
                hold on
            else

                scatter(figHandle, Signal(:,1), Signal(:,2), 36, 'filled');
                hold on
                if ~isempty(Noise)
                    scatter(figHandle, Noise(:,1), Noise(:,2), 20, 'k', 'filled', 'MarkerEdgeColor','w');
                    hold on

                end
            end

            if ~strcmp(obj.Chains(iChains).name, "V")
                gm = obj.Chains(iChains).gmmModel.Model;

                theta = linspace(0,2*pi,100);
                conf = chi2inv(0.95,2);                % scale for 95% ellipse
                for j = 1:gm.NumComponents
                    Sigma = gm.Sigma(:,:,j);
                    [V,D] = eig(Sigma);
                    a = sqrt(conf*diag(D));            % axis lengths
                    ellipse = (V*(a.*[cos(theta); sin(theta)])) + gm.mu(j,:)';
                    plot(figHandle, ellipse(1,:), ellipse(2,:), 'LineWidth',1.5);
                end
            end
            
            title([obj.Chains(iChains).name, ' Measure Data'])

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
                % Extract data from chain                
                txData.Chains(iChains).dataType  = obj.Chains(iChains).name;
                txData.Chains(iChains).data      = obj.Chains(iChains).gmmModel.filteredData;
                txData.Chains(iChains).Means     = obj.Chains(iChains).gmmModel.Model.mu;
                % txData.Chains(iChains).Variances = diag(obj.Chains(iChains).gmmModel.Model.Sigma(:,:,1));
                txData.Chains(iChains).Variances = obj.Chains(iChains).gmmModel.Model.Sigma;
                txData.Chains(iChains).prob      = obj.Chains(iChains).gmmModel.Model.ComponentProportion;
                txData.Chains(iChains).Corrupt   = obj.Chains(iChains).gmmModel.badDataIdx;
                
            end
        end
    
    end

end