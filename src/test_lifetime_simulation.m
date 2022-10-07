function test_lifetime_simulation()
    % Run lifetime simulation with the discrete SEI growth model
    %
    % 2022/08/03
    % - just demonstrate how we can use the discrete update version of the
    %   SEI growth model to handle changing use cases.


    % Initialize the system
    p = initialize_cell_props();

    % Define end condition
    dQa_vec = 0;
    dRa_vec = 0;
    dQa1_vec = 0;
    dRa1_vec = 0;
    dQa2_vec = 0;
    dRa2_vec = 0;

    % Initialize loop variables
    preva.dQ = 0;
    preva.dR = 0;
    preva1.dQ = 0;
    preva1.dR = 0;
    preva2.dQ = 0;
    preva2.dR = 0;

    MAX_CYCLE_NUMBER = 1000;
    cyc_num_vec = 0;
    cyc_num = 1;

    za0 = 0;
    zb0 = 0;
    alpha = 1.2;
    U0 = 3.0;
    Vmin = 3.0;
    Vmax = 4.2;

    f_ocv = @(z) U0 + alpha * z;

    hours_to_charge = 3;
    current_target = -p.Qb / (hours_to_charge * 3600);
    t = linspace(0, 10*3600, 1.0e6)';
    I_chg = +current_target*ones(size(t)); % applied current (A) 
    I_dch = -current_target*ones(size(t));
    I_current_cutoff = current_target/20;

    % Loop over cycles
    while cyc_num < MAX_CYCLE_NUMBER
    
        fprintf('Simulating cycle %g ...\n', cyc_num)

        % Update capacity loss and resistance growth increment for this cycle
        delta_time = hours_to_charge * 3600;

        curra = update_states_simple(preva, delta_time, cyc_num);
        curra1 = update_states_simple(preva1, delta_time, 1);
        curra2 = update_states_simple(preva2, delta_time, 900);

        dQa_vec = [dQa_vec ; curra.dQ];
        dRa_vec = [dRa_vec ; curra.dR];
        dQa1_vec = [dQa1_vec ; curra1.dQ];
        dRa1_vec = [dRa1_vec ; curra1.dR];
        dQa2_vec = [dQa2_vec ; curra2.dQ];
        dRa2_vec = [dRa2_vec ; curra2.dR];
        preva = curra;
        preva1 = curra1;
        preva2 = curra2;

        cyc_num_vec = [cyc_num_vec; cyc_num];
        cyc_num = cyc_num + 1;

    end

    set_default_plot_settings()

    figure('Position', [10 10 700 600].*1.5);
    subplot(211)
    line(cyc_num_vec, (p.Qa - dQa1_vec)./3600, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Usage A');
    line(cyc_num_vec, (p.Qa - dQa2_vec)./3600, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Usage B');
    line(cyc_num_vec, (p.Qa - dQa_vec)./3600, 'Color', 'k', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 1.5, 'DisplayName', 'Usage A $\rightarrow$ B');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('Capacity (Ah)', 'Interpreter', 'Latex')
    ylim([0 5])
    lh = legend('show');

    subplot(212)
    line(cyc_num_vec, (p.Ra + dRa1_vec).*1000, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Usage A');
    line(cyc_num_vec, (p.Ra + dRa2_vec).*1000, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Usage B');
    line(cyc_num_vec, (p.Ra + dRa_vec).*1000, 'Color', 'k', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 1.5, 'DisplayName', 'Usage A $\rightarrow$ B');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('Resistance (m$\Omega$)', 'Interpreter', 'Latex')
    lh = legend('show', 'Interpreter', 'Latex');

end

function curr = update_states_simple(prev, delta_time, cyc_num)
    % Q-R state update function 
    %
    % Parameters
    % ---------
    % prev (struct)
    %    previous state (includes dQ and dR)
    % delta_time (float)
    %    time for the cycle
    % Iss (float)
    %    steady-state current value in Amperes
    %
    % Outputs
    % ---------
    % curr
    %    updated state

    n = 0.5;

    % Actually we're going to change the k after some number of cycles,
    % just to demonstrate the idea of dynamically updating k values based
    % on changing use cases
    if cyc_num < 500
        kc = 1;
        kr = 0.00001;
    else 
        kc = 2;
        kr = 0.00002;
    end
    
    % Discrete state update equation for SEI growth
    curr.dQ = (kc^(1/n) * delta_time + prev.dQ^(1/n))^n;
    curr.dR = (kr^(1/n) * delta_time + prev.dR^(1/n))^n;

end