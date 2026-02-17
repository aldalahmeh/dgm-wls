classdef GridEnvironment < handle
    %GridEnvironment Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Grid
        nMonteCarlo
        ThermalNoiseParam
        ImpNoiseParam        
        TrueGridState
    end

    methods
        function obj = GridEnvironment(Grid, ThermalNoiseParam, ImpNoiseParam, nMonteCarlo)
            %GridEnvironment Construct an instance of this class
            %   Detailed explanation goes here
            obj.Grid = Grid;
            obj.TrueGridState = Grid.TrueGridState;
            obj.ThermalNoiseParam = ThermalNoiseParam;
            obj.ImpNoiseParam = ImpNoiseParam;
            obj.nMonteCarlo = nMonteCarlo;
        end

        function obj = genSensorMeas(obj)
            %genSensorMeas Summary of this method goes here
            %   Detailed explanation goes here
            
            nBus = numel(obj.Grid.Bus);
            measSelector = obj.Grid.sensorTable;

            for iBus = 1:nBus
                % Add noise for each sensor measurement taken from iBus
                noisyThermalMeas = obj.addThermalNoise(iBus, measSelector);
                impulsiveNoise = obj.addImpulsiveNoise(iBus, measSelector);
                
                obj.Grid.Bus(iBus).meas = noisyThermalMeas + impulsiveNoise;
            end
        end

        function noisyMeas = addThermalNoise(obj, iBus, measSelector)                                

            % avalDataType = obj.getAvalDataType(measSelector);
            % basVal = obj.Grid.Bus(iBus).trueVal{:,avalDataType};
            % noisStd = basVal .* obj.ThermalNoiseParam{:,avalDataType};

            % Get sensor data types for iBus
            dataTypeInd = ( measSelector.Bus == iBus);
            % dataType = measSelector.DataType(dataTypeInd);   
            datatype = measSelector.DataType(dataTypeInd);
            branchStr = "-"+ string(measSelector.Branch(dataTypeInd)); 
            branchStr(ismissing(branchStr)) = "";

            dataType = datatype + branchStr; 

            % val = table2array(measSelector);
            % Get base values and noise std for such sensor data types
            % basVal = obj.Grid.Bus(iBus).trueVal.(val);
            % noisStd = basVal .* obj.ThermalNoiseParam.(val);
            basVal = obj.Grid.Bus(iBus).trueVal(:,dataType);
            basVal = table2array(basVal);

            noisStd = basVal .* obj.ThermalNoiseParam(:,extractBefore(dataType,3));
            noisStd = table2array(noisStd);

            noisyMeasIntr = basVal + noisStd .* randn(obj.nMonteCarlo, length(basVal)) ;          
            noisyMeas = array2table(noisyMeasIntr, 'VariableNames', dataType);

           
        end

        function implMeas = addImpulsiveNoise(obj, iBus, measSelector)

            pImp = obj.ImpNoiseParam.prob;

            % val = table2array(measSelector);

            % avalDataType = obj.getAvalDataType(measSelector);
            
            % basVal = obj.Grid.Bus(iBus).trueVal{:,avalDataType};
            % basVal = obj.Grid.Bus(iBus).trueVal.(val);

            % Get sensor data types for iBus
            dataTypeInd = ( measSelector.Bus == iBus);
            % dataType = measSelector.DataType(dataTypeInd);
            % dataType = measSelector.DataType(dataTypeInd) + "-" ...
            %     + string(measSelector.Branch(dataTypeInd)); 
            datatype = measSelector.DataType(dataTypeInd);
            branchStr = "-"+ string(measSelector.Branch(dataTypeInd)); 
            branchStr(ismissing(branchStr)) = "";

            dataType = datatype + branchStr; 

            basVal = obj.Grid.Bus(iBus).trueVal(:,dataType);
            basVal = table2array(basVal);

            
            ImpNoisIndicator = rand(obj.nMonteCarlo,length(basVal)) < pImp;
            noisStd = basVal .* obj.ImpNoiseParam.level;
            implMeasIntr = ImpNoisIndicator .* (basVal + noisStd .* randn(obj.nMonteCarlo, length(basVal)));
         
            implMeas = array2table(implMeasIntr, 'VariableNames', dataType);
        end

        function avalDataType = getAvalDataType(obj, measSelector)

            fn = fieldnames(measSelector);             % cell array of field names
            vals = cellfun(@(f) measSelector.(f), fn); % logical array (same order as fn)
            avalDataType = fn(vals);

        end
    end
end