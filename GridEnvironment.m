classdef GridEnvironment < handle
    %GRIDENVIRONMENT - Environment wrapper for grid simulations
    %   OBJ = GRIDENVIRONMENT(GRID, THERMALNOISEPARAM, IMPNOISEPARAM, NSIM)
    %   constructs an environment that holds grid configuration and noise
    %   parameters used across NSIM simulation runs.
    %
    %   Properties:
    %       Grid               - grid configuration object
    %       nSim               - number of simulation realizations
    %       ThermalNoiseParam  - parameters for thermal noise model
    %       GridSensorNoiseStd - per-sensor noise std storage (cell)
    %       ImpNoiseParam      - impulse noise parameters
    %       TrueGridState      - reference true grid state
    % 
    %   See also HANDLE
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
            % Initialize core environment fields from inputs
            obj.Grid = Grid;
            obj.TrueGridState = Grid.TrueGridState;
            obj.ThermalNoiseParam = ThermalNoiseParam;
            obj.GridSensorNoiseStd = {};
            obj.ImpNoiseParam = ImpNoiseParam;
            obj.nSim = nSim;
        end
        
        %% Getters
        function impNoiseProp = getImpNoiseProp(obj)
          % GETIMPNOISEPROP - Return impulse noise probability from parameters
            impNoiseProp = obj.ImpNoiseParam.prob;
        end
        function nSim = getNumSim(obj)
          % GETNUMSIM - Return configured number of simulations
            nSim = obj.nSim;
        end
        function impNoiseMultiplier = getImpNoiseMultiplier(obj)
          % GETIMPNOISEMULTIPLIER - Return impulse noise amplitude multiplier
            impNoiseMultiplier = obj.ImpNoiseParam.multiplier;
        end
   
        function sensorNoiseStd = getSensorNoiseStd(obj, dataType)
          % GETSENSORNOISESTD - Return thermal noise std; optional per-dataType slice
          %
          % Input arguments:
          % obj      - object containing ThermalNoiseParam
          % dataType - (optional) string whose prefix selects columns
            if nargin == 1
                sensorNoiseStd = obj.ThermalNoiseParam;
            elseif nargin == 2
            sensorNoiseStd = obj.ThermalNoiseParam(:,extractBefore(dataType,3));
            end
        end
        function busImpNoisIndicator = getBusImpNoisInd(obj, iBus)
          % GETBUSIMPNOISIND - Retrieve impedance noise indicator for a bus
            busImpNoisIndicator = obj.Grid.Bus(iBus).impNoisIndicator;
        end

        function gridImpNoisIndicator = getGridImpNoisInd(obj, testDataRange)
          % GETGRIDIMPNOISIND - Extract per-bus impedance/noise indicators for the requested range
          %
          % Input arguments:
          % obj           - object containing Grid with Bus entries
          % testDataRange - row indices (time/samples) to extract
          %
          % Output arguments:
          % gridImpNoisIndicator - cell array: each cell is subset for one bus
            nBus = numel(obj.Grid.Bus);
            for iBus = 1:nBus
              % Get full indicator matrix for bus iBus, then subset rows
                intrImpInd = getBusImpNoisInd(obj, iBus);
                gridImpNoisIndicator{iBus} = intrImpInd(testDataRange,:);
            end
        end
        
        function gridSensorNoiseStd = getGridSensorNoiseStd(obj)
          % GETGRIDSENSORNOISESTD - Assemble per-bus sensor noise std tables
          %
          % Input arguments:
          % obj - object providing Grid and sensor metadata
          %
          % Output arguments:
          % gridSensorNoiseStd - cell array of per-bus tables of sensor std
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
         
        end
        %%
        function obj = genSensorMeas(obj)
          % GENSENSORMEAS - Generate sensor measurements for each bus
          %
          % Input arguments:
          % obj - object containing Grid and sensor configuration
          %
          % Output arguments:
          % obj - object updated with noisy measurements and noise metadata
            %genSensorMeas Summary of this method goes here
            %   Detailed explanation goes here
            
            nBus = numel(obj.Grid.Bus);
            measSelector = obj.Grid.sensorTable;
            for iBus = 1:nBus
                % Add noise for each sensor measurement taken from iBus
                [noisyThermalMeas, thermalStd] = obj.addThermalNoise(iBus, measSelector);
                [impulsiveNoise, implInd] = obj.addImpulsiveNoise(iBus, measSelector, thermalStd);
                
                obj.Grid.Bus(iBus).meas = noisyThermalMeas + impulsiveNoise;
                
                % Store wich measurement got corrupted
                obj.Grid.Bus(iBus).impNoisIndicator = implInd;
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
            noisMeasStd = table2array(obj.getSensorNoiseStd(dataType));
           
           
            % Compute sensor noise floor
            % The max value is 3*sigma, so divide by 3
            noiseFloorStd = 1e-4;

            % Compare with noise floor, choose the max
            noisStd = max(noisMeasStd, noiseFloorStd);

            % Create noisy measurements
            noisyMeasIntr = measNominalValPU + noisStd .* randn(nSim, length(measNominalVal)) ;   

            % Convert to tables
            noisyMeas = array2table(noisyMeasIntr, 'VariableNames', dataType);
            noisStdTable = array2table(noisStd, 'VariableNames', dataType);
           
        end

        function [implMeas, implInd] = addImpulsiveNoise(obj, iBus, measSelector, thermalStd)
          % ADDIMPULSIVENOISE - Add impulsive noise samples to selected measurements
          %
          % Input arguments:
          % obj          - object providing simulation and grid data
          % iBus         - bus index for which measurements are considered
          % measSelector - logical/index selector for measurements at the bus
          % thermalStd   - thermal noise std table for measurement types
          %
          % Output arguments:
          % implMeas - table of impulsive noise samples (nSim x nMeas)
          % implInd  - binary indicator matrix for impulsive occurrences
            % Number of simulation data points
            nSim = obj.getNumSim();          
            % Get imulsive noise level
            impNoiseMultiplier = obj.getImpNoiseMultiplier();
            % Get measured data types at bus
            dataType = obj.getDataTypes(iBus, measSelector);
            % Get nominal values of measurements
            measNominalVal = obj.Grid.getTrueVal(iBus, dataType);
                              
            % Randomly set impulsive noise indicator at sensor-couples
            [impNoisIndicator, implInd] = obj.setImpNoiseIndSensors(dataType);
            % DEBUGGING
            % DEBUG = true;
            DEBUG = false;
            if DEBUG 
                % Force impulses at specific rows for debugging temporal behavior
                impNoisIndicator(2001,:) = ones(1,length(dataType));
                impNoisIndicator(2002:end,:) = zeros(1999,length(dataType));            
            end
             
            % Compute standard deviation, multiplier x thermal std 
            noisStd = table2array(thermalStd .* impNoiseMultiplier);
            % Impulsive measurement array
            % Multiply indicator by std and Gaussian noise to generate impulses
            implMeasIntr = impNoisIndicator .* noisStd .* randn( nSim, length(measNominalVal) );
         
            % Impulsive measurement table
            implMeas = array2table(implMeasIntr, 'VariableNames', dataType);
           
        end

        %% Helper functions

        function perUnitmeas = normalizePerUnit(obj, meas, chainDataType)
          % NORMALIZEPERUNIT - Convert measurements to grid per-unit where appropriate
          %
          % Input arguments:
          % obj - object providing Grid.getBaseVal()
          % meas - numeric array of measurements
          % chainDataType - string/array indicating measurement types (e.g., "VM")
          
            % Grid base value
            basVal = obj.Grid.getBaseVal();
           % Normalize measurements to per-unit
           perUnitmeas = meas / basVal;
           vmInd = strcmp(chainDataType, "VM");
           perUnitmeas(vmInd) = meas(vmInd);
           
        end
        function sensorCouples = mapMeastoSensCouple(~, dataType)
          % MAPMEASTOSENSCOUPLE - Group measurement types into sensor couples
          %
          % Input arguments:
          % ~        - unused object handle
          % dataType - string array of measurement type identifiers
          %
          % Output arguments:
          % sensorCouples - string array of grouped sensor identifiers
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
                    % Pair active and reactive power measurements into one sensor
                    sensorCouples(iCell) = {"PD-QD"};
                    idataType = idataType + 2;
                elseif strcmp(dataType(idataType), "PG") && strcmp(dataType(idataType+1), "QG")
                    % Pair generator active/reactive measurements into one sensor
                    sensorCouples(iCell) = {"PG-QG"};
                    idataType = idataType + 2;
                elseif contains(dataType(idataType), "PF") && contains(dataType(idataType+1), "QF")
                    % Pair feeder flow with corresponding reactive measurement, preserve suffix
                    sensorCouples(iCell) = {"PF-" + dataType(idataType+1)};
                    idataType = idataType + 2;
                end
                iCell = iCell + 1;
            end
            % Convert to string array
            sensorCouples = string([sensorCouples{:}]);
        end
        function [impNoisIndicator, impSensorNoisIndicatorTable] ...
                = setImpNoiseIndSensors(obj, dataType)
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
            impSensorNoisIndicatorTable = array2table(ImpSensorNoisIndicatorIntr,'VariableNames',sensorCouples);
            
            % Copy inticator to measurements in pairs
            impNoisIndicator = zeros(nSim, length(dataType));

            % Number of available data types
            nSensorCouples = numel(sensorCouples);
               
            iCol = 1;           % column index counter
        
            % Go through the columns of ImpNoisIndicator
            for iSensorCouples = 1:nSensorCouples
                % Set the column as the corresponding column in ImpSensorNoisIndicator
                 if strcmp(sensorCouples(iSensorCouples), "VM") 
                    impNoisIndicator(:,iCol) = impSensorNoisIndicatorTable.(sensorCouples(iSensorCouples));
                    iCol = iCol + 1;
                 
                 % Set the column as the corresponding column in ImpSensorNoisIndicator 
                 % Copy the same ImpSensorNoisIndicator column to the
                 % adjacent column in the case of PD-QD, PG-QG and PF-QF
                 else
                    impNoisIndicator(:,iCol) = impSensorNoisIndicatorTable.(sensorCouples(iSensorCouples));
                    impNoisIndicator(:,iCol+1) = impSensorNoisIndicatorTable.(sensorCouples(iSensorCouples));
                    iCol = iCol + 2;

                 end                 
            end                       

        end
        function DataTypes = getDataTypes(~, iBus, measSelector)
          % GETDATATYPES - Build combined data type strings for a bus
          % Input arguments:
          % iBus - bus index to filter measurements
          % measSelector - struct with fields Bus, DataType, Branch
            % Get sensor data types for iBus
            dataTypeInd = ( measSelector.Bus == iBus);
            datatype = measSelector.DataType(dataTypeInd);
            branchStr = "-"+ string(measSelector.Branch(dataTypeInd)); 
            branchStr(ismissing(branchStr)) = "";
            DataTypes = datatype + branchStr; 
        end
        function avalDataType = getAvalDataType(~, measSelector)
          % GETAVALDATATYPE - Return field names whose values are true
          % measSelector is a struct of logical flags; returns names set true.
            fn = fieldnames(measSelector);             % cell array of field names
            vals = cellfun(@(f) measSelector.(f), fn); % logical array (same order as fn)
            avalDataType = fn(vals);
        end
    end
end