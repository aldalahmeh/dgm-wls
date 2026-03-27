classdef MonteCarloPerformance
    % MONTECARLOPERFORMANCE Evaluates, processes, and plots Monte Carlo simulation results.
    %   This class provides a comprehensive suite of tools for computing Mean 
    %   Absolute Error (MAE), convergence statistics, and detection metrics (F1-score) 
    %   from power grid state estimation data. It also includes methods to generate 
    %   publication-ready figures and LaTeX tables.
    
    properties
        resultsPath % Directory path where simulation results (.mat files) are saved and loaded
    end
    
    methods
        %% Constructor
        function obj = MonteCarloPerformance(filePath)
            % MONTECARLOPERFORMANCE Constructs the performance evaluation object.
            obj.resultsPath = filePath;
        end
        
        %% Core Computation Methods
        function maeMetrics = computeMAEfromFile(obj, data)
            % COMPUTEMAEFROMFILE Extracts results from a loaded struct and computes MAE.
            results  = data.Algorithms.MAEMetrics;
            trueVals = data.TrueGridState;
            maeMetrics = obj.computeMAE(results, trueVals);
        end
        
        function maeMetrics = computeMAE(obj, results, trueVals)
            % COMPUTEMAE Computes the Mean Absolute Error for Voltage Magnitude and Angle.
            % Creates matrices with VM/VA estimates across all Monte Carlo runs.
            vmEstMat = obj.createMat(results, "VM");
            vaEstMat = obj.createMat(results, "VA");
            
            % Calculate MAE by averaging across the Monte Carlo iterations (dimension 2)
            maeMetrics.VM = mean(abs(vmEstMat - trueVals.VM), 2);
            maeMetrics.VA = mean(abs(vaEstMat - trueVals.VA), 2);
        end
        
        function convMetrics = computeConvStats(~, convData)
            % COMPUTECONVSTATS Calculates convergence rate and average iterations.
            nSim = length(convData);
            nonConvNum = sum(isnan(convData));
            
            % Metrics calculation
            convMetrics.convRate = (nSim - nonConvNum) / nSim;
            convMetrics.aveIter = mean(convData(~isnan(convData)));
        end
        
        function f1Score = computeF1Score(obj, trueImp, estImp)
            % COMPUTEF1SCORE Calculates the F1 detection score for impulsive noise across the grid.
            nBus = numel(trueImp);
            TP = 0; % True Positive
            TN = 0; % True Negative
            FP = 0; % False Positive
            FN = 0; % False Negative
            
            % Aggregate detection statistics across all buses
            for iBus = 1:nBus
                [busTP, busTN, busFP, busFN] = obj.cmptDetStat(trueImp{iBus}, estImp{iBus});
                TP = TP + busTP;
                TN = TN + busTN;
                FP = FP + busFP;
                FN = FN + busFN;
            end
            
            % Calculate Precision safely to avoid division by zero
            Precision = TP / max((TP + FP), realmin);
            % Calculate Recall safely to avoid division by zero
            Recall = TP / max((TP + FN), realmin);
            
            % Calculate the final F1 Score (Harmonic mean)
            if (Precision + Recall) == 0
                f1Score = 0;
            else
                f1Score = 2 * (Precision * Recall) / (Precision + Recall);
            end
        end
        
        function [busTP, busTN, busFP, busFN] = cmptDetStat(~, busTrueImp, busEstImp)
            % CMPTDETSTAT Computes True Positives, False Positives, etc., for a specific bus.
            measTypes = busTrueImp.Properties.VariableNames;
            nMeasType = numel(measTypes);
            
            for iMeasType = 1:nMeasType
                trueImpInd = busTrueImp.(char(measTypes(iMeasType)));
                estImpInd = busEstImp.(char(measTypes(iMeasType)));
                
                % TP: Ground truth and estimate identify an anomaly
                busTP = sum(trueImpInd & estImpInd);
                % TN: Ground truth and estimate identify clean data
                busTN = sum(~trueImpInd & ~estImpInd);
                % FP: Ground truth is clean, but estimate flags an anomaly (False Alarm)
                busFP = sum(~trueImpInd & estImpInd);
                % FN: Ground truth is anomaly, but estimate flags as clean (Missed Detection)
                busFN = sum(trueImpInd & ~estImpInd);
            end
        end
        
        %% Plotting Methods
        function plotCompareMAE(~, varargin)
            % PLOTCOMPAREMAE Plots a side-by-side grouped bar chart comparison for multiple methods.
            % Usage: obj.plotCompareMAE(mae1, mae2, ..., title1, title2, ...)
            
            nInputs = length(varargin);
            if mod(nInputs, 2) ~= 0
                error('Inputs must be paired: N metric structs followed by N titles.');
            end
            
            nMethods = nInputs / 2;
            maeMetrics = varargin(1:nMethods);
            figTitles = varargin(nMethods+1:end);
            
            % Extract data into matrices
            nBus = length(maeMetrics{1}.VM);
            vmData = zeros(nBus, nMethods);
            vaData = zeros(nBus, nMethods);
            
            for i = 1:nMethods
                vmData(:, i) = maeMetrics{i}.VM(:);
                vaData(:, i) = maeMetrics{i}.VA(:);
            end
            
            colors = lines(nMethods);
            
            % --- Figure for VM Comparison ---
            figure;
            b_vm = bar(vmData, 'grouped');
            for i = 1:nMethods
                b_vm(i).FaceColor = colors(i, :);
            end
            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title('Mean Absolute Error for VM Estimates');
            legend(figTitles{:}, 'Location', 'north', 'FontSize', 12);
            grid on;
            
            % --- Figure for VA Comparison ---
            figure;
            b_va = bar(vaData, 'grouped');
            for i = 1:nMethods
                b_va(i).FaceColor = colors(i, :);
            end
            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title('Mean Absolute Error for VA Estimates');
            legend(figTitles{:}, 'Location', 'north', 'FontSize', 12);
            grid on;
        end
        
        function plotMAE(~, maeMetrics, figTitle)
            % PLOTMAE Plots individual MAE bar charts for a single algorithm.
            figure;
            bar(maeMetrics.VM, 'FaceColor', 'b');
            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title(['Mean Absolute Error for VM Estimates - ' figTitle]);
            legend('VM');
            grid on;
            
            figure;
            bar(maeMetrics.VA, 'FaceColor', 'r');
            xlabel('Bus Index');
            ylabel('Mean Absolute Error');
            title(['Mean Absolute Error for VA Estimates- ' figTitle]);
            legend('VA');
            grid on;
        end
        
        function plotConvStats(~, varargin)
            % PLOTCONVSTATS Plots the convergence rate and average iterations comparison.
            nInputs = length(varargin);
            if mod(nInputs, 2) ~= 0
                error('Inputs must be paired: N metric structs followed by N titles.');
            end
            
            nMethods = nInputs / 2;
            convMetrics = varargin(1:nMethods);
            figTitles = varargin(nMethods+1:end);
            
            successRates = zeros(1, nMethods);
            avgIters = zeros(1, nMethods);
            
            for i = 1:nMethods
                successRates(i) = convMetrics{i}.convRate * 100;
                avgIters(i)     = convMetrics{i}.aveIter;
            end
            
            colors = lines(nMethods);
            
            % --- Figure 1: Convergence Rate ---
            figure; hold on;
            for i = 1:nMethods
                b1 = bar(i, successRates(i));
                b1.FaceColor = colors(i, :);
            end
            hold off;
            ylabel('Convergence Rate (%)');
            title('Convergence Reliability');
            set(gca, 'XTick', 1:nMethods, 'XTickLabel', figTitles);
            ylim([0 105]); 
            grid on;
            
            % --- Figure 2: Average Iterations ---
            figure; hold on;
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
        
        function plotLoadedDataMAE(obj, Data)
            % PLOTLOADEDDATAMAE Computes and plots MAE directly from a loaded .mat struct.
            % Enforces IEEE-compliant formatting including strict Y-axis ticks for VM plots.
            trueVals = Data.TrueGridState;
            nMethods = length(Data.Algorithms);
            
            maeMetricsCell = cell(1, nMethods);
            titlesCell = cell(1, nMethods);
            
            for i = 1:nMethods
                if isfield(Data.Algorithms(i), 'Name')
                    titlesCell{i} = Data.Algorithms(i).Name;
                else
                    titlesCell{i} = sprintf('Method %d', i);
                end
                
                if isfield(Data.Algorithms(i), 'StateEst')
                    estResults = Data.Algorithms(i).StateEst;
                elseif isfield(Data.Algorithms(i), 'results')
                    estResults = Data.Algorithms(i).results;
                else
                    error('Could not find StateEst or results field in the Algorithms struct.');
                end
                
                maeMetricsCell{i} = obj.computeMAE(estResults, trueVals);
            end
            
            plotArgs = [maeMetricsCell, titlesCell];
            obj.plotCompareMAE(plotArgs{:});
            
            % Adjust Y-Axis Ticks Strictly for VM across all open figures
            allAxes = findall(0, 'Type', 'axes', 'Tag', '');
            for k = 1:length(allAxes)
                ax = allAxes(k);
                ax.FontSize = 10;
                
                if ischar(ax.Title.String)
                    titleStr = ax.Title.String;
                elseif iscell(ax.Title.String) && ~isempty(ax.Title.String)
                    titleStr = ax.Title.String{1};
                else
                    titleStr = '';
                end
                
                if contains(titleStr, 'VM')
                    yLimits = ylim(ax);
                    tickStep = 0.001; 
                    ax.YTick = 0 : tickStep : yLimits(2);
                    ax.YGrid = 'on';
                else
                    ax.YTickMode = 'auto';
                    ax.YGrid = 'on';
                end
            end
            
            % Enforce global legend styling
            allLegends = findall(0, 'Type', 'legend');
            if ~isempty(allLegends)
                set(allLegends, 'FontSize', 12, 'Location', 'north');
            end
            
            fprintf('Successfully plotted and formatted MAE comparison for %d algorithms.\n', nMethods);
        end
        
        function plotAvgIterTrend(obj, fileNames, probabilities)
            % PLOTAVGITERTREND Plots the average iterations trend across varying noise probabilities.
            nFiles = length(fileNames);
            if nFiles ~= length(probabilities)
                error('Number of files must match the number of probabilities.');
            end
            
            firstData = obj.load(fileNames{1});
            nAlgs = length(firstData.Algorithms);
            algNames = cell(1, nAlgs);
            
            for j = 1:nAlgs
                if isfield(firstData.Algorithms(j), 'Name')
                    algNames{j} = firstData.Algorithms(j).Name;
                else
                    algNames{j} = sprintf('Method %d', j);
                end
            end
            
            % Override names for consistency
            if nAlgs >= 3
                algNames{1} = 'WLS';
                algNames{2} = 'DGM-WLS';
                algNames{3} = 'IRWLS';
            end
            
            avgItersMat = zeros(nFiles, nAlgs);
            for i = 1:nFiles
                Data = obj.load(fileNames{i});
                for j = 1:nAlgs
                    if isfield(Data.Algorithms(j), 'ConvMetrics') && isfield(Data.Algorithms(j).ConvMetrics, 'aveIter')
                        avgItersMat(i, j) = Data.Algorithms(j).ConvMetrics.aveIter;
                    elseif isfield(Data.Algorithms(j), 'ConvData')
                        convStats = obj.computeConvStats(Data.Algorithms(j).ConvData);
                        avgItersMat(i, j) = convStats.aveIter;
                    else
                        error('Could not find convergence data for algorithm %s in file %s', algNames{j}, fileNames{i});
                    end
                end
            end
            
            figure;
            b = bar(avgItersMat, 'grouped');
            colors = lines(nAlgs);
            for j = 1:nAlgs
                b(j).FaceColor = colors(j, :);
            end
            
            probLabels = cell(1, nFiles);
            for i = 1:nFiles
                probLabels{i} = sprintf('p = %g', probabilities(i));
            end
            set(gca, 'XTick', 1:nFiles, 'XTickLabel', probLabels);
            xlabel('Impulsive Error Probability');
            ylabel('Average Number of Iterations');
            title('Computational Effort vs. Noise Probability');
            legend(algNames{:}, 'Location', 'northeast');
            grid on;
            
            fprintf('Successfully plotted Average Iterations bar chart across %d scenarios.\n', nFiles);
        end
        
        function plotAvgMAETrend(obj, fileNames, probabilities)
            % PLOTAVGMAETREND Plots the average grid-wide MAE trend for both VM and VA.
            nFiles = length(fileNames);
            if nFiles ~= length(probabilities)
                error('Number of files must match the number of probabilities.');
            end
            
            firstData = obj.load(fileNames{1});
            nAlgs = length(firstData.Algorithms);
            algNames = cell(1, nAlgs);
            
            for j = 1:nAlgs
                if isfield(firstData.Algorithms(j), 'Name')
                    algNames{j} = firstData.Algorithms(j).Name;
                else
                    algNames{j} = sprintf('Method %d', j);
                end
            end
            
            if nAlgs >= 3
                algNames{1} = 'WLS';
                algNames{2} = 'DGM-WLS';
                algNames{3} = 'IRWLS';
            end
            
            avgMaeVmMat = zeros(nFiles, nAlgs);
            avgMaeVaMat = zeros(nFiles, nAlgs);
            
            for i = 1:nFiles
                Data = obj.load(fileNames{i});
                trueVals = Data.TrueGridState;
                
                for j = 1:nAlgs
                    if isfield(Data.Algorithms(j), 'StateEst')
                        estResults = Data.Algorithms(j).StateEst;
                    elseif isfield(Data.Algorithms(j), 'results')
                        estResults = Data.Algorithms(j).results;
                    else
                        error('Could not find StateEst or results field in the Algorithms struct.');
                    end
                    
                    maeMetrics = obj.computeMAE(estResults, trueVals);
                    avgMaeVmMat(i, j) = mean(maeMetrics.VM);
                    avgMaeVaMat(i, j) = mean(maeMetrics.VA);
                end
            end
            
            probLabels = cell(1, nFiles);
            for i = 1:nFiles
                probLabels{i} = sprintf('p = %g', probabilities(i));
            end
            colors = lines(nAlgs);
            
            % --- VM MAE Trend ---
            figure;
            b_vm = bar(avgMaeVmMat, 'grouped');
            for j = 1:nAlgs
                b_vm(j).FaceColor = colors(j, :);
            end
            set(gca, 'XTick', 1:nFiles, 'XTickLabel', probLabels);
            xlabel('Impulsive Error Probability');
            ylabel('Average MAE (VM)');
            title('Average VM MAE vs. Noise Probability');
            legend(algNames{:}, 'Location', 'northwest');
            grid on;
            
            % --- VA MAE Trend ---
            figure;
            b_va = bar(avgMaeVaMat, 'grouped');
            for j = 1:nAlgs
                b_va(j).FaceColor = colors(j, :);
            end
            set(gca, 'XTick', 1:nFiles, 'XTickLabel', probLabels);
            xlabel('Impulsive Error Probability');
            ylabel('Average MAE (VA)');
            title('Average VA MAE vs. Noise Probability');
            legend(algNames{:}, 'Location', 'northwest');
            grid on;
            
            fprintf('Successfully plotted Average VM and VA MAE trends across %d scenarios.\n', nFiles);
        end
        
        %% LaTeX Formatting
        function generateF1LatexTable(obj, matFileName, texFileName)
            % GENERATEF1LATEXTABLE Reads tuning results and outputs a publication-ready LaTeX table.
            data = obj.load(matFileName);
            minpts = data.minptsValues;
            epsVals = data.epsilonValues;
            f1Mat = data.f1Score;
            maxF1 = max(f1Mat(:));
            
            texFileName = string(texFileName);
            if ~endsWith(texFileName, ".tex", 'IgnoreCase', true)
                texFileName = texFileName + ".tex";
            end
            
            texFilePath = fullfile(obj.resultsPath, texFileName);
            fid = fopen(texFilePath, 'w');
            if fid == -1
                error('Could not open file %s for writing.', texFilePath);
            end
            
            nCols = length(epsVals);
            colFormat = repmat('c', 1, nCols);
            
            fprintf(fid, '\\begin{table}[h]\n');
            fprintf(fid, '\\centering\n');
            fprintf(fid, '\\caption{DBSCAN Hyperparameter Tuning: $F_1$ Scores}\n');
            fprintf(fid, '\\label{tab:dbscan_tuning}\n');
            fprintf(fid, '\\begin{tabular}{@{}l%s@{}}\n', colFormat);
            fprintf(fid, '\\toprule\n');
            fprintf(fid, 'MinPts / $\\epsilon$');
            
            for j = 1:nCols
                fprintf(fid, ' & \\textbf{%.2f}', epsVals(j));
            end
            fprintf(fid, ' \\\\\n\\midrule\n');
            
            for i = 1:length(minpts)
                fprintf(fid, '\\textbf{%d}', minpts(i));
                for j = 1:nCols
                    val = f1Mat(i, j);
                    if abs(val - maxF1) < 1e-6
                        fprintf(fid, ' & \\textbf{%.7f}', val);
                    else
                        fprintf(fid, ' & %.7f', val);
                    end
                end
                fprintf(fid, ' \\\\\n');
            end
            
            fprintf(fid, '\\bottomrule\n');
            fprintf(fid, '\\end{tabular}\n');
            fprintf(fid, '\\end{table}\n');
            fclose(fid);
            
            fprintf('\n--- Generated LaTeX Table ---\n\n');
            type(texFilePath);
            fprintf('\n-----------------------------\n');
            fprintf('Table successfully saved to: %s\n', texFilePath);
        end
        
        %% File I/O & Export
        function save(obj, fileName, algSelect, Algorithms, varargin)
            % SAVE Stores simulation configuration, results, and custom variables to a .mat file.
            fileName = string(fileName);
            if ~endsWith(fileName, ".mat", 'IgnoreCase', true)
                fileName = fileName + ".mat";
            end
            fullPath = fullfile(obj.resultsPath, fileName);
            
            save(fullPath, "algSelect", "Algorithms");
            savedVars = {};
            
            if ~isempty(varargin)
                if ischar(varargin{1}) || isstring(varargin{1})
                    savedVars = cellstr(string(varargin));
                    args = sprintf(',''%s''', savedVars{:});
                    evalin('caller', sprintf('save(''%s''%s, ''-append'')', fullPath, args));
                else
                    savedVars = cell(1, numel(varargin));
                    for k = 1:numel(varargin)
                        vName = inputname(k + 4);
                        if isempty(vName)
                            error('Arguments must be named variables in the workspace.');
                        end
                        savedVars{k} = vName;
                        eval(sprintf('%s = varargin{k};', vName));
                    end
                    save(fullPath, savedVars{:}, '-append');
                end
            end
            
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
        
        function saveF1Score(~, fileName, minptsValues, epsilonValues, f1Score)
            % SAVEF1SCORE Exports hyperparameter tuning variables to a MAT file.
            save(fileName, 'minptsValues', 'epsilonValues', 'f1Score');
        end
        
        function Results = load(obj, fileName)
            % LOAD Retrieves results from a designated .mat file.
            fullPath = fullfile(obj.resultsPath, fileName);
            Results = load(fullPath);
            fprintf('Data loaded from file: %s\n', fullPath);
        end
        
        function saveFigures(~, figsPath, simID, fileExt)
            % SAVEFIGURES Formats open figures to IEEE standards and exports them.
            if ~exist(figsPath, 'dir')
                mkdir(figsPath);
            end
            fileExt = strrep(fileExt, '.', '');
            figs = findall(0, 'Type', 'figure');
            
            if isempty(figs)
                fprintf('No open figures found to save.\n');
                return;
            end
            
            fprintf('\n--- Saving Figures ---\n');
            for i = 1:length(figs)
                fig = figs(i);
                set(fig, 'Color', 'w');
                ax = findall(fig, 'Type', 'axes');
                
                if isempty(ax)
                    continue;
                end
                
                set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', 'k');
                set(findall(fig, 'Type', 'text'), 'Color', 'k');
                set(findall(fig, 'Type', 'legend'), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
                
                titleStr = string(ax(1).Title.String);
                
                if contains(titleStr, "Magnitude", "IgnoreCase", true) || contains(titleStr, "VM", "IgnoreCase", true)
                    prefix = "mae_vm";
                elseif contains(titleStr, "Angle", "IgnoreCase", true) || contains(titleStr, "VA", "IgnoreCase", true)
                    prefix = "mae_va";
                elseif contains(titleStr, "Reliability", "IgnoreCase", true) || contains(titleStr, "Rate", "IgnoreCase", true)
                    prefix = "conv_rate";
                elseif contains(titleStr, "Effort", "IgnoreCase", true) || contains(titleStr, "Iteration", "IgnoreCase", true)
                    prefix = "conv_iter";
                else
                    prefix = sprintf("figure_%d", fig.Number);
                end
                
                fileName = sprintf('%s%s.%s', prefix, simID, fileExt);
                fullPath = fullfile(figsPath, fileName);
                
                for j = 1:length(ax)
                    ax(j).Title.String = '';
                    if isprop(ax(j), 'Subtitle')
                        ax(j).Subtitle.String = '';
                    end
                end
                
                delete(findall(fig, 'Tag', 'sgtitle'));
                delete(findall(fig, 'Tag', 'suptitle'));
                
                if strcmpi(fileExt, 'eps')
                    exportgraphics(fig, fullPath, 'ContentType', 'vector', 'BackgroundColor', 'w');
                else
                    exportgraphics(fig, fullPath, 'Resolution', 300, 'BackgroundColor', 'w');
                end
                fprintf('Saved: %s\n', fullPath);
            end
            fprintf('----------------------\n');
        end
        
        function saveTargetFigure(~, targetFig, figsPath, fileName)
            % SAVETARGETFIGURE Invisibly formats and saves a single specific target figure.
            if ~exist(figsPath, 'dir')
                mkdir(figsPath);
            end
            if isnumeric(targetFig)
                targetFig = figure(targetFig);
            end
            
            fprintf('\n--- Saving Target Figure ---\n');
            hiddenFig = copyobj(targetFig, 0);
            set(hiddenFig, 'Visible', 'off', 'Color', 'w');
            ax = findall(hiddenFig, 'Type', 'axes');
            
            if ~isempty(ax)
                set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', 'GridColor', 'k');
                for j = 1:length(ax)
                    ax(j).Title.Color = 'k';
                    if isprop(ax(j), 'Subtitle') && ~isempty(ax(j).Subtitle)
                        ax(j).Subtitle.Color = 'k';
                    end
                end
            end
            
            set(findall(hiddenFig, 'Tag', 'sgtitle'), 'Color', 'k');
            legs = findall(hiddenFig, 'Type', 'legend');
            if ~isempty(legs)
                set(legs, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
            end
            
            allObjs = findall(hiddenFig);
            colorProps = {'Color', 'MarkerEdgeColor', 'MarkerFaceColor', 'TextColor', 'EdgeColor'};
            for i = 1:length(allObjs)
                obj = allObjs(i);
                if strcmpi(obj.Type, 'figure') || strcmpi(obj.Type, 'axes') || strcmpi(obj.Type, 'legend')
                    continue;
                end
                for p = 1:length(colorProps)
                    prop = colorProps{p};
                    if isprop(obj, prop)
                        try
                            val = obj.(prop);
                            isWhite = isequal(val, [1 1 1]) || ...
                                (ischar(val) && (strcmpi(val, 'w') || strcmpi(val, 'white'))) || ...
                                (isstring(val) && (strcmpi(val, "w") || strcmpi(val, "white")));
                            if isWhite
                                obj.(prop) = 'k';
                            end
                        catch
                            % Ignore read-only properties
                        end
                    end
                end
            end
            
            fullPath = fullfile(figsPath, fileName);
            [~, ~, fileExt] = fileparts(fileName);
            if strcmpi(fileExt, '.eps') || strcmpi(fileExt, '.pdf')
                exportgraphics(hiddenFig, fullPath, 'ContentType', 'vector', 'BackgroundColor', 'w');
            else
                exportgraphics(hiddenFig, fullPath, 'Resolution', 300, 'BackgroundColor', 'w');
            end
            delete(hiddenFig);
            fprintf('Successfully saved: %s\n', fullPath);
            fprintf('----------------------------\n');
        end
        
        %% Utility functions
        function vxEstMat = createMat(~, dataCellArray, colName)
            % CREATEMAT Aggregates Monte Carlo estimate arrays into a single unified matrix.
            nMonteCarlo = numel(dataCellArray);
            nBus = height(dataCellArray{1});
            emptyInd = [];
            vxEstMat = zeros(nBus, nMonteCarlo);
            
            for i = 1:nMonteCarlo
                if isempty(dataCellArray{i}) || any(isnan(dataCellArray{i}.VM)) || any(isnan(dataCellArray{i}.VA))
                    emptyInd = [emptyInd, i];
                    continue;
                else
                    vxEstMat(:, i) = dataCellArray{i}.(colName);
                end
            end
            
            vxEstMat(:, emptyInd) = []; 
        end
    end
end