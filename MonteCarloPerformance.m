classdef MonteCarloPerformance
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        resultsPath
    end

    methods
        function obj = MonteCarloPerformance(filePath)
            %UNTITLED Construct an instance of this class
            obj.resultsPath = filePath;

        end

        function maeMetrics = computeMAEfromFile(obj, data)

            results  = data.results;
            trueVals = data.trueVals;

            maeMetrics = obj.computeMAE(results, trueVals);

        end

        function maeMetrics = computeMAE(obj, results, trueVals)
            % COMPUTENAE Summary of this method goes here
            %   Detailed explanation goes here



            % Create a matrix with VM estimates
            vmEstMat = obj.createMat(results, "VM");
            vaEstMat = obj.createMat(results, "VA");

            % Calculate Mean Absolute Error (MAE)
            maeMetrics.VM = mean(abs(vmEstMat - trueVals.VM), 2);
            maeMetrics.VA = mean(abs(vaEstMat - trueVals.VA), 2);


        end

        function convMetrics = computeConvStats(~, convData)

            % Number of simulation iterations
            nSim = length(convData);

            % Number of non-convergence iterations
            nonConvNum = sum(isnan(convData));

            % Convergence rate
            convMetrics.convRate = (nSim - nonConvNum)/nSim;

            % Average number of convergence iterations
            convMetrics.aveIter = mean(convData(~isnan(convData)));

        end
        %% Plotting
        function plotCompareMAE(~, varargin)
            % PLOTCOMPAREMAE Plots a side-by-side comparison for an arbitrary number of methods
            % Usage: obj.plotCompareMAE(mae1, mae2, ..., title1, title2, ...)

            % 1. Determine the number of inputs and split them
            nInputs = length(varargin);
            if mod(nInputs, 2) ~= 0
                error('Inputs must be paired: N metric structs followed by N titles.');
            end

            nMethods = nInputs / 2;
            maeMetrics = varargin(1:nMethods);
            figTitles = varargin(nMethods+1:end);

            % 2. Dynamically extract data into N x M matrices
            % Get number of buses from the first struct to preallocate
            nBus = length(maeMetrics{1}.VM);
            vmData = zeros(nBus, nMethods);
            vaData = zeros(nBus, nMethods);

            for i = 1:nMethods
                vmData(:, i) = maeMetrics{i}.VM(:);
                vaData(:, i) = maeMetrics{i}.VA(:);
            end

            % 3. Generate a dynamic color palette
            % 'lines' creates MATLAB's default highly distinguishable colors
            colors = lines(nMethods);

            % --- Figure for VM Comparison ---
            figure;
            b_vm = bar(vmData, 'grouped');

            % Apply dynamic styling
            for i = 1:nMethods
                b_vm(i).FaceColor = colors(i, :);
            end

            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title('Mean Absolute Error for VM Estimates');
            legend(figTitles{:}, 'Location', 'best');
            grid on;

            % --- Figure for VA Comparison ---
            figure;
            b_va = bar(vaData, 'grouped');

            % Apply dynamic styling
            for i = 1:nMethods
                b_va(i).FaceColor = colors(i, :);
            end

            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title('Mean Absolute Error for VA Estimates');
            legend(figTitles{:}, 'Location', 'best');
            grid on;
        end

        function plotMAE(~, maeMetrics, figTitle)

            % Figure for VM
            figure;
            bar(maeMetrics.VM, 'FaceColor', 'b');
            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title(['Mean Absolute Error for VM Estimates - ' figTitle]);
            legend('VM');
            grid on;

            % Figure for VA
            figure;
            bar(maeMetrics.VA, 'FaceColor', 'r');
            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title(['Mean Absolute Error for VA Estimates- ' figTitle]);
            legend('VA');
            grid on;


        end

        function plotConvStats(~, varargin)

            % 1. Determine the number of inputs and split them
            nInputs = length(varargin);
            if mod(nInputs, 2) ~= 0
                error('Inputs must be paired: N metric structs followed by N titles.');
            end

            nMethods = nInputs / 2;
            convMetrics = varargin(1:nMethods);
            figTitles = varargin(nMethods+1:end);

            % Preallocate arrays for the metrics
            successRates = zeros(1, nMethods);
            avgIters = zeros(1, nMethods);

            % 2. Extract the pre-calculated metrics for each method
            for i = 1:nMethods
                successRates(i) = convMetrics{i}.convRate * 100;
                avgIters(i)     = convMetrics{i}.aveIter;
            end

            % Generate dynamic color palette
            colors = lines(nMethods);

            % --- Figure 1: Convergence Rate ---
            figure;
            hold on;
            for i = 1:nMethods
                b1 = bar(i, successRates(i));
                b1.FaceColor = colors(i, :);
            end
            hold off;

            ylabel('Convergence Rate (%)');
            title('Convergence Reliability');
            set(gca, 'XTick', 1:nMethods, 'XTickLabel', figTitles);
            ylim([0 105]); % Cap at 100% with slight headroom for visibility
            grid on;

            % --- Figure 2: Average Iterations ---
            figure;
            hold on;
            for i = 1:nMethods
                b2 = bar(i, avgIters(i));
                b2.FaceColor = colors(i, :);
            end
            hold off;

            ylabel('Average Number of Iterations');
            title('Computational Effort (Successful Runs Only)');
            set(gca, 'XTick', 1:nMethods, 'XTickLabel', figTitles);
            grid on;

        end

        %% File I/O
        function save(obj, fileName, algSelect, Algorithms, varargin)
            % SAVE Stores the simulation configuration, algorithms struct,
            % and any additional requested variables to a .mat file.

            % 1. Format the filename and path
            fileName = string(fileName);
            if ~endsWith(fileName, ".mat", 'IgnoreCase', true)
                fileName = fileName + ".mat";
            end
            fullPath = fullfile(obj.resultsPath, fileName);

            % 2. Save the Core Metadata & Results
            % Because Algorithms contains the StateEst, ConvMetrics, and Options,
            % saving it alongside algSelect guarantees perfect reproducibility.
            save(fullPath, "algSelect", "Algorithms");

            % Initialize list of additionally saved variables for reporting
            savedVars = {};

            % 3. Handle any additional variables (varargin)
            if ~isempty(varargin)
                % Check if the user passed variable names as strings (e.g., "trueVals")
                if ischar(varargin{1}) || isstring(varargin{1})
                    savedVars = cellstr(string(varargin));
                    % Construct comma-separated string of variable names for evalin
                    args = sprintf(',''%s''', savedVars{:});
                    % Append to the existing file
                    evalin('caller', sprintf('save(''%s''%s, ''-append'')', fullPath, args));

                    % Check if the user passed the actual variables (e.g., trueVals)
                else
                    savedVars = cell(1, numel(varargin));
                    for k = 1:numel(varargin)
                        % inputname offsets: 1=obj, 2=fileName, 3=algSelect, 4=Algorithms
                        vName = inputname(k + 4);
                        if isempty(vName)
                            error('Arguments must be named variables in the workspace.');
                        end
                        savedVars{k} = vName;

                        % Assign value to the local workspace with its true name
                        eval(sprintf('%s = varargin{k};', vName));
                    end
                    % Append only the dynamically named variables
                    save(fullPath, savedVars{:}, '-append');
                end
            end

            % 4. Generate the Save Summary Report
            % Extract the names of the algorithms that were executed and saved
            if ~isempty(Algorithms)
                savedAlgs = strjoin([Algorithms.Name], ', ');
            else
                savedAlgs = 'None';
            end

            fprintf('\n--- Save Summary ---\n');
            fprintf('File: %s\n', fullPath);
            fprintf('Core Data Saved: algSelect, Algorithms\n');
            fprintf('Algorithms Included: %s\n', savedAlgs);
            if ~isempty(savedVars)
                fprintf('Additional Variables: %s\n', strjoin(savedVars, ', '));
            end
            fprintf('--------------------\n');
        end

        function Results = load(obj, fileName)

            fullPath = fullfile(obj.resultsPath, fileName);
            Results = load(fullPath);

            fprintf('Data loaded from file: %s\n', fullPath);
        end

        function saveFigures(~, figsPath, simID, fileExt)
            % SAVEFIGURES Automatically detects open figures, categorizes them,
            % and saves them in publication-ready quality.

            % Ensure the directory exists
            if ~exist(figsPath, 'dir')
                mkdir(figsPath);
            end

            % Clean up the extension if the user passed it with a dot
            fileExt = strrep(fileExt, '.', '');

            % Get all currently open figures
            figs = findall(0, 'Type', 'figure');

            if isempty(figs)
                fprintf('No open figures found to save.\n');
                return;
            end

            fprintf('\n--- Saving Figures ---\n');
            for i = 1:length(figs)
                fig = figs(i);

                % Find the axes to read the title
                ax = findobj(fig, 'Type', 'axes');
                if isempty(ax)
                    continue; % Skip empty figures
                end

                % Grab the title string of the first axis
                titleStr = string(ax(1).Title.String);

                % Smart Categorization based on Title content
                if contains(titleStr, "Magnitude", "IgnoreCase", true) || contains(titleStr, "VM", "IgnoreCase", true)
                    prefix = "mae_vm";
                elseif contains(titleStr, "Angle", "IgnoreCase", true) || contains(titleStr, "VA", "IgnoreCase", true)
                    prefix = "mae_va";
                elseif contains(titleStr, "Reliability", "IgnoreCase", true) || contains(titleStr, "Rate", "IgnoreCase", true)
                    prefix = "conv_rate";
                elseif contains(titleStr, "Effort", "IgnoreCase", true) || contains(titleStr, "Iteration", "IgnoreCase", true)
                    prefix = "conv_iter";
                else
                    % Fallback for unknown figures
                    prefix = sprintf("figure_%d", fig.Number);
                end

                % Construct full filename
                fileName = sprintf('%s%s.%s', prefix, simID, fileExt);
                fullPath = fullfile(figsPath, fileName);

                % Export settings based on file type
                if strcmpi(fileExt, 'eps')
                    % EPS requires vector format for crisp IEEE scaling
                    exportgraphics(fig, fullPath, 'ContentType', 'vector', 'BackgroundColor', 'none');
                else
                    % PNG/JPEG get 300 DPI high-resolution
                    exportgraphics(fig, fullPath, 'Resolution', 300);
                end

                fprintf('Saved: %s\n', fullPath);
            end
            fprintf('----------------------\n');
        end

        %% Utility functions
        function vxEstMat = createMat(~, dataCellArray, colName)


            nMonteCarlo = numel(dataCellArray);
            nBus = height(dataCellArray{1});

            emptyInd = [];
            % Initialize the matrix based on the number of Monte Carlo points
            vxEstMat = zeros(nBus, nMonteCarlo);
            for i = 1:nMonteCarlo
                if isempty(dataCellArray{i}) || any(isnan(dataCellArray{i}.VM)) || any(isnan(dataCellArray{i}.VA))
                    % Save empty data indicies
                    emptyInd = [emptyInd, i];
                    continue
                else
                    vxEstMat(:, i) = dataCellArray{i}.(colName);
                end
            end

            % Remove empty (zero) columns
            vxEstMat(:, emptyInd) = []; % Remove empty columns from the matrix

        end
    end
end