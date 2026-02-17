function busSensorSelector = create_nonzero_sensor_selection(casename)
    % create_nonzero_sensor_selection
    % Creates a sensor selection struct based on NON-ZERO values in the case file.
    % This ensures we only simulate sensors where physical power actually exists.
    
    % 1. Define Constants
    define_constants;
    mpc = loadcase(casename);

    % 2. Initialize
    nb = size(mpc.bus, 1);
    template = struct('VM',false, 'PD',false, 'QD',false, ...
                      'PF',false, 'QF',false, 'PG',false, 'QG',false);
    busSensorSelector = repmat(template, nb, 1);
    
    % 3. Pre-process Generation Data
    % Generators might not be ordered by bus, so we sum them per bus first.
    % (Some buses might have multiple gens, though rare in case14)
    P_Gen_Total = sparse(mpc.gen(:, GEN_BUS), 1, mpc.gen(:, PG), nb, 1);
    Q_Gen_Total = sparse(mpc.gen(:, GEN_BUS), 1, mpc.gen(:, QG), nb, 1);
    
    % 4. Pre-process Branch Data (Sending Ends)
    % Find which buses are the "From" side of any branch
    sending_end_buses = unique(mpc.branch(:, F_BUS));
    
    % 5. Loop through each bus
    for i = 1:nb
        
        % --- A. Voltage Sensors (VM) ---
        % Logic: Measure V only at Voltage-Controlled Buses (PV or Slack)
        % Bus Types: 3=Ref, 2=PV.
        if mpc.bus(i, BUS_TYPE) == 3 || mpc.bus(i, BUS_TYPE) == 2
            busSensorSelector(i).VM = true;
        end
        
        % --- B. Load Sensors (PD, QD) ---
        % Logic: Only place a sensor if there is actual Load defined (> 1e-4 MW)
        % This filters out Bus 1, Bus 7, Bus 8.
        if abs(mpc.bus(i, PD)) > 1e-4 || abs(mpc.bus(i, QD)) > 1e-4
            busSensorSelector(i).PD = true;
            busSensorSelector(i).QD = true;
        end
        
        % --- C. Generation Sensors (PG, QG) ---
        % Logic: Only place a sensor if there is active/reactive Generation
        if abs(P_Gen_Total(i)) > 1e-4 && abs(Q_Gen_Total(i)) > 1e-4
            busSensorSelector(i).PG = true;
            busSensorSelector(i).QG = true;
        end
        
        % --- D. Flow Sensors (PF, QF) ---
        % Logic: Measure Sending End flows if lines start here.
        if ismember(i, sending_end_buses)
            busSensorSelector(i).PF = true;
            busSensorSelector(i).QF = true;
        end
    end
end