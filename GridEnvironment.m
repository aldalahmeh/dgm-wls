classdef GridEnvironment < handle
    %GridEnvironment Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Grid
        nSim
        ThermalNoiseParam
        GridSensorNoiseStd
        ImpNoiseParam        
        TrueGridState
    end

    methods
        %% Constructor
        function obj = GridEnvironment(Grid, ThermalNoiseParam, ImpNoiseParam, nSim)

            obj.Grid = Grid;
            obj.TrueGridState = Grid.TrueGridState;
            obj.ThermalNoiseParam = ThermalNoiseParam;
            obj.GridSensorNoiseStd = {};
            obj.ImpNoiseParam = ImpNoiseParam;
            obj.nSim = nSim;
        end
        
        %% Getters
        function impNoiseProp = getImpNoiseProp(obj)
            
            impNoiseProp = obj.ImpNoiseParam.prob;
        end

        function nSim = getNumSim(obj)
            
            nSim = obj.nSim;
        end

        function impNoiseMultiplier = getImpNoiseMultiplier(obj)

            impNoiseMultiplier = obj.ImpNoiseParam.multiplier;
        end

        % function sensorError = getSensorError(obj, dataType)
        % 
        %     sensorError = obj.ThermalNoiseParam(:,extractBefore(dataType,3));
        % 
        % end

        function sensorNoiseStd = getSensorNoiseStd(obj, dataType)

            if nargin == 1

                sensorNoiseStd = obj.ThermalNoiseParam;

            elseif nargin == 2

            sensorNoiseStd = obj.ThermalNoiseParam(:,extractBefore(dataType,3));
            end

        end

        % !!
        function gridSensorNoiseStd = getGridSensorNoiseStd(obj)

            % baseVal = obj.Grid.getBaseVal();
            sensorErrorTable = obj.getSensorNoiseStd();

            nBus = obj.Grid.getNumBus();
            measSelector = obj.Grid.getMeasSel();

            for iBus = 1:nBus 

                % Get measurements data types of the sensors at iBus
                dataType = obj.getDataTypes(iBus, measSelector);

                % Get corresponding data-type error from sensor table error
                T = sensorErrorTable(:,extractBefore(dataType,3));

                % Correct the PF/QF names in the table header.
                % !!
                colNames = string(T.Properties.VariableNames);
                mask = startsWith(colNames, "PF") | startsWith(colNames, "QF");
                T.Properties.VariableNames(:,mask) = dataType(mask)';

                gridSensorNoiseStd{iBus} = T;

            end

            % sensorNoiseStd = {};
            % nCells = numel(sensorNoiseStdCellAry);
            % for iCell = 1:nCells
            %     T = sensorNoiseStdCellAry{iCell};
            %     % vars = T.Properties.VariableNames;        % names
            %     % varsToScale = vars(~strcmp(vars,'VM'));   % exclude 'VM'
            %     % if ~isempty(varsToScale)
            %     %     T{:,varsToScale} = T{:,varsToScale} / baseVal;  % scale numeric data
            %     % end
            %     sensorNoiseStd{iCell} = T;
            % end

        end
        %%
        function obj = genSensorMeas(obj)
            %genSensorMeas Summary of this method goes here
            %   Detailed explanation goes here
            
            nBus = numel(obj.Grid.Bus);
            measSelector = obj.Grid.sensorTable;

            for iBus = 1:nBus
                % Add noise for each sensor measurement taken from iBus
                [noisyThermalMeas, thermalStd] = obj.addThermalNoise(iBus, measSelector);
                impulsiveNoise = obj.addImpulsiveNoise(iBus, measSelector, thermalStd);
                
                obj.Grid.Bus(iBus).meas = noisyThermalMeas + impulsiveNoise;
                % !!
                obj.GridSensorNoiseStd(iBus) = {thermalStd};
            end
        end
     
        %%
        function [noisyMeas, noisStdTable] = addThermalNoise(obj, iBus, measSelector)                                

                 
            % Number of simulation data points
            nSim = obj.getNumSim();
           
            % Get measurements data types of the sensors at iBus
            dataType = obj.getDataTypes(iBus, measSelector);

            % Get nominal measurements values of the sensors at iBus
            measNominalVal = obj.Grid.getTrueVal(iBus, dataType);

            % Normalize to pu 
            measNominalValPU = obj.normalizePerUnit(measNominalVal, dataType);

            % Compute measurement standard deviation                        
            % sensorError = table2array(obj.getSensorError(dataType));
            noisMeasStd = table2array(obj.getSensorNoiseStd(dataType));
           
            % noisMeasStd = (abs(measNominalVal) .* sensorError)/3;
            % noisMeasStd = (abs(measNominalVal) .* sensorError);
            
            % Get measurements max. values from matpower mpc
            % fullScale = obj.Grid.getFullScale(iBus, dataType);

            % Compute sensor noise floor
            % The max value is 3*sigma, so divide by 3
            % noiseFloorStd = ((sensorError/10) .* fullScale)/3;
            % noiseFloorStd = ((sensorError/10) .* fullScale);
            noiseFloorStd = 1e-4;

            % Compare with noise floor, choose the max
            noisStd = max(noisMeasStd, noiseFloorStd);

            % Create noisy measurements
            noisyMeasIntr = measNominalValPU + noisStd .* randn(nSim, length(measNominalVal)) ;   

            % Convert to tables
            noisyMeas = array2table(noisyMeasIntr, 'VariableNames', dataType);
            noisStdTable = array2table(noisStd, 'VariableNames', dataType);
           
        end

        function implMeas = addImpulsiveNoise(obj, iBus, measSelector, thermalStd)

            
            % Number of simulation data points
            nSim = obj.getNumSim();          

            % Get imulsive noise level
            impNoiseMultiplier = obj.getImpNoiseMultiplier();

            % Get measured data types at bus
            dataType = obj.getDataTypes(iBus, measSelector);

            % Get nominal values of measurements
            measNominalVal = obj.Grid.getTrueVal(iBus, dataType);
                              
            % Randomly set impulsive noise indicator at sensor-couples
            impNoisIndicator = obj.setImpNoiseIndSensors(dataType);

            % DEBUGGING
            % DEBUG = true;
            DEBUG = false;

            if DEBUG 
                % Set bus measurements to be corrupt for all buses only for
                % the first iteration and not the rest of iteration
                impNoisIndicator(2001,:) = ones(1,length(dataType));
                impNoisIndicator(2002:end,:) = zeros(1999,length(dataType));            
            end
             
            % Compute standard deviation, multiplier x thermal std 
            noisStd = table2array(thermalStd .* impNoiseMultiplier);

            % Impulsive measurement array
            implMeasIntr = impNoisIndicator .* noisStd .* randn( nSim, length(measNominalVal) );
         
            % Impulsive measurement table
            implMeas = array2table(implMeasIntr, 'VariableNames', dataType);
        end

        %% Helper functions

        function perUnitmeas = normalizePerUnit(obj, meas, chainDataType)
          
            % Grid base value
            basVal = obj.Grid.getBaseVal();

           % Normalize measurements to per-unit
           perUnitmeas = meas / basVal;

           vmInd = strcmp(chainDataType, "VM");

           perUnitmeas(vmInd) = meas(vmInd);

            % if strcmp(chainDataType, "VM")
            %     perUnitmeas = meas;
            % else
            %     perUnitmeas = meas / basVal;
            % end

        end

        function sensorCouples = mapMeastoSensCouple(~, dataType)

            % Initialize cell
            sensorCouples = {};
            % Number of available data types
            ndataType = numel(dataType);
               
            idataType = 1; % data type counter
            iCell = 1;     % cell index counter
            while idataType <= ndataType
                if strcmp(dataType(idataType), "VM")
                    sensorCouples(iCell) = {"VM"};
                    idataType = idataType + 1;
                    

                elseif strcmp(dataType(idataType), "PD") && strcmp(dataType(idataType+1), "QD")
                    sensorCouples(iCell) = {"PD-QD"};
                    idataType = idataType + 2;

                elseif strcmp(dataType(idataType), "PG") && strcmp(dataType(idataType+1), "QG")
                    sensorCouples(iCell) = {"PG-QG"};
                    idataType = idataType + 2;

                elseif contains(dataType(idataType), "PF") && contains(dataType(idataType+1), "QF")
                    sensorCouples(iCell) = {"PF-" + dataType(idataType+1)};
                    idataType = idataType + 2;

                end
                iCell = iCell + 1;
            end

            % Convert to string array
            sensorCouples = string([sensorCouples{:}]);


        end
        function ImpNoisIndicator = setImpNoiseIndSensors(obj, dataType)
            %IMPULSEINDICATOR returns a matrix of logical values indicating
            % the occurance of an impulsive noise at a specific bus having
            % sensor measurements determined by input argument DATATYPE as 
            % the columns of IMPULSEINDICATOR whereas the rows are
            % represent the Monte Carlo iteration.

            % Number of simulation data points
            nSim = obj.getNumSim();

            % Get imulsive noise probability                       
            pImp = obj.getImpNoiseProp();

            % Map measurements to sensor-couples
            sensorCouples = obj.mapMeastoSensCouple(dataType);
            
            % Randomly set sensors indicator
            ImpSensorNoisIndicatorIntr = rand(nSim, length(sensorCouples)) < pImp;

            % Create table for ease of access
            ImpSensorNoisIndicator = array2table(ImpSensorNoisIndicatorIntr,'VariableNames',sensorCouples);
            
            % Copy inticator to measurements in pairs
            ImpNoisIndicator = zeros(nSim, length(dataType));

            % Number of available data types
            nSensorCouples = numel(sensorCouples);
               
            % iSensorCouples = 1; % data type counter
            iCol = 1;           % column index counter
        
            % Go through the columns of ImpNoisIndicator
            for iSensorCouples = 1:nSensorCouples
                % Set the column as the corresponding column in ImpSensorNoisIndicator
                 if strcmp(sensorCouples(iSensorCouples), "VM") 
                    ImpNoisIndicator(:,iCol) = ImpSensorNoisIndicator.(sensorCouples(iSensorCouples));
                    iCol = iCol + 1;
                 
                 % Set the column as the corresponding column in ImpSensorNoisIndicator 
                 % Copy the same ImpSensorNoisIndicator column to the
                 % adjacent column in the case of PD-QD, PG-QG and PF-QF
                 else
                    ImpNoisIndicator(:,iCol) = ImpSensorNoisIndicator.(sensorCouples(iSensorCouples));
                    ImpNoisIndicator(:,iCol+1) = ImpSensorNoisIndicator.(sensorCouples(iSensorCouples));
                    iCol = iCol + 2;

                 end

                 % iSensorCouples = iSensorCouples + 1;
                 
            end                       

        end
        function DataTypes = getDataTypes(~, iBus, measSelector)

            % Get sensor data types for iBus
            dataTypeInd = ( measSelector.Bus == iBus);

            datatype = measSelector.DataType(dataTypeInd);
            branchStr = "-"+ string(measSelector.Branch(dataTypeInd)); 
            branchStr(ismissing(branchStr)) = "";

            DataTypes = datatype + branchStr; 
        end
        function avalDataType = getAvalDataType(~, measSelector)

            fn = fieldnames(measSelector);             % cell array of field names
            vals = cellfun(@(f) measSelector.(f), fn); % logical array (same order as fn)
            avalDataType = fn(vals);

        end
    end
end