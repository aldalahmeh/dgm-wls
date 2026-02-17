classdef CommsChannel
    %COMMSCHANNEL Summary of this class goes here
    %   Detailed explanation goes here

    properties
        SNR_dB
    end

    methods
        function obj = CommsChannel(SNR_dB)
            % COMMSCHANNEL Construct an instance of this class
            %   Detailed explanation goes here
            obj.SNR_dB = SNR_dB;
        end

        function rxData = addWhiteNoise(obj, txData)
            % ADDWHITENOISE Adds AWGN based on a specified SNR (Signal-to-Noise Ratio)
            
            rxData = txData;
            nChains = numel(txData.Chains);
            
            % Assume you have defined obj.SNR_dB (e.g., 30 dB)
            snr_linear = 10^(obj.SNR_dB / 10);
            
            for iChains = 1:nChains
                rxData.Chains(iChains).dataType = txData.Chains(iChains).dataType;
                data = txData.Chains(iChains).data;
                [nData, nDataType] = size(data);
                
                % 1. Calculate Average Signal Power for this specific measurement type
                % Using root mean square (RMS) power of the dataset
                signalPower = mean(data.^2, 'all'); 
                
                % Fallback for zero-power signals (e.g., Zero Injection Buses)
                if signalPower < 1e-6
                    signalPower = 1e-4; % Assume a small base power
                end
                
                % 2. Calculate the fixed Communication Variance based on SNR
                commVariance = signalPower / snr_linear;
                commStdDev = sqrt(commVariance);
                
                % 3. Generate independent AWGN
                noise = commStdDev .* randn(nData, nDataType);
                
                % 4. Store the true variance for the WLS W-Matrix
                % Notice this is now a constant value for the whole chain, 
                % correctly representing a stable communication link.
                rxData.Chains(iChains).commVar = repmat(commVariance, nData, nDataType);
                
                % 5. Add noise to the clean data
                rxData.Chains(iChains).data = data + noise;
            end            
        end
        function rxData = addWhiteNoiseFlawed(obj, txData)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            rxData = txData;
            nChains = numel(txData.Chains);
            
            % Define a tiny noise floor to prevent Zero-Variance matrices
            noiseFloor = 1e-6; 
            
            for iChains = 1:nChains
                rxData.Chains(iChains).dataType = txData.Chains(iChains).dataType;
                [nData, nDataType] = size(txData.Chains(iChains).data);
                
                % 1. Calculate Standard Deviation (Use abs() for magnitude)
                % 2. Add noiseFloor to prevent 0-variance on zero-injection/zero-flow data
                basVal = (obj.noiseStd .* abs(txData.Chains(iChains).data)) + noiseFloor;
                
                % Generate scaled Gaussian noise
                noise = basVal .* randn(nData, nDataType);
                
                % Store the true variance (sigma^2)
                rxData.Chains(iChains).commVar = basVal.^2;
                
                % Add noise to the clean data
                rxData.Chains(iChains).data = txData.Chains(iChains).data + noise;
            end        
        end
    end
end