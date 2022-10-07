function test_lifetime_simulation_with_imbalance()
    % Run lifetime simulation for two parallel-connected batteries
    %
    % 2022/07/29
    % - Use the imbalance dynamics equations to simulate intra-cycle dynamics
    % - At the end of each cycle, compute aggregate metrics and us it to
    %   update the capacities and resistances.
    % - Repeat.


    % Initialize the system
    p = initialize_cell_props();

    % Define end condition
    dQa_vec = 0;
    dQb_vec = 0;
    dRa_vec = 0;
    dRb_vec = 0;
    dQa_control_vec = 0;
    dQb_control_vec = 0;
    dRa_control_vec = 0;
    dRb_control_vec = 0;
    Issa_vec = nan;
    Issb_vec = nan;
    zssa_vec = nan;
    zssb_vec = nan;

    % Initialize loop variables
    preva.dQ = 0;
    prevb.dQ = 0;
    preva.dR = 0;
    prevb.dR = 0;

    % Loop variables for the control cells (no current sharing)
    preva_control.dQ = 0;
    prevb_control.dQ = 0;
    preva_control.dR = 0;
    prevb_control.dR = 0;

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

        % Solve for the intra-cycle dynamics
        %
        % curr: struct holding current (Ra, Rb, Qa, Rb)
        % prev: struct holding previous (Ra, Rb, Qa, Rb)
        Qa = p.Qa - preva.dQ;
        Qb = p.Qb - prevb.dQ;
        Ra = p.Ra + preva.dR;
        Rb = p.Rb + prevb.dR;

        % Update capacity loss and resistance growth increment for this cycle
        delta_time = hours_to_charge * 3600;

        [Issa, Issb, zssa, zssb] = update_cycle_metrics(Qa, Qb, Ra, Rb, alpha, I_chg(1), za0, zb0, delta_time);
        
        curra = update_states_iss(preva, delta_time, Issa);
        currb = update_states_iss(prevb, delta_time, Issb);

        curra_control = update_states_iss(preva_control, delta_time, I_chg(1)/2);
        currb_control = update_states_iss(prevb_control, delta_time, I_chg(1)/2);

        dQa_vec = [dQa_vec ; curra.dQ];
        dQb_vec = [dQb_vec ; currb.dQ];
        dRa_vec = [dRa_vec ; curra.dR];
        dRb_vec = [dRb_vec ; currb.dR];

        dQa_control_vec = [dQa_control_vec ; curra_control.dQ];
        dQb_control_vec = [dQb_control_vec ; currb_control.dQ];
        dRa_control_vec = [dRa_control_vec ; curra_control.dR];
        dRb_control_vec = [dRb_control_vec ; currb_control.dR];  

        Issa_vec = [Issa_vec ; Issa];
        Issb_vec = [Issb_vec ; Issb];
        zssa_vec = [zssa_vec ; zssa];
        zssb_vec = [zssb_vec ; zssb];

        preva = curra;
        prevb = currb;
        preva_control = curra_control;
        prevb_control = currb_control;

        cyc_num_vec = [cyc_num_vec; cyc_num];
        cyc_num = cyc_num + 1;

    end

    set_default_plot_settings()

    figure('Position', [10 10 700 600].*1.5);
    subplot(221)
    line(cyc_num_vec, (p.Qa - dQa_vec)./3600, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec, (p.Qa - dQa_control_vec)./3600, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', ':');
    line(cyc_num_vec, (p.Qb - dQb_vec)./3600, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec, (p.Qb - dQb_control_vec)./3600, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', ':');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('Capacity (Ah)', 'Interpreter', 'Latex')
    ylim([0 5])
    lh = legend('show');

    subplot(222)
    line(cyc_num_vec, (p.Ra + dRa_vec).*1000, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec, (p.Ra + dRa_control_vec).*1000, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', ':');
    line(cyc_num_vec, (p.Rb + dRb_vec).*1000, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec, (p.Rb + dRb_control_vec).*1000, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', ':');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('Resistance (m$\Omega$)', 'Interpreter', 'Latex')
    lh = legend('show');

    subplot(223)
    line(cyc_num_vec, Issa_vec, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec, Issb_vec, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('$I_{ss}$ (A)', 'Interpreter', 'Latex')
    lh = legend('show');

    subplot(224)
    line(cyc_num_vec, zssa_vec, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec, zssb_vec, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('$z_{ss}$', 'Interpreter', 'Latex')
    lh = legend('show');

end

function [Issa, Issb, zssa, zssb] = update_cycle_metrics(Qa, Qb, Ra, Rb, alpha, I, za0, zb0, dt)
    % Compute steady-state values

    Q = Qa + Qb;
    R = Ra + Rb;

    kappa = 1/alpha * (Ra*Qa - Rb*Qb) / Q;

    Issa = 1/R * ((Rb + alpha*kappa))*I;
    Issb = 1/R * ((Ra - alpha*kappa))*I;
    
    zssa = 1/Q * (Qa*za0 + Qb*zb0) + 1/Q * (+Qb*kappa - dt) * I;
    zssb = 1/Q * (Qa*za0 + Qb*zb0) + 1/Q * (-Qa*kappa - dt) * I;

    % For future consideration, we can also run a full discrete-time
    % simulation here:

%         out = run_discrete_time_simulation_complete(I_chg, I_dch, ...
%                cv_cutoff_current_amps, Qa, Qb, Ra, Rb, za0, zb0, f_ocv, ...
%                Vmin, Vmax);

%         out = solve_z_dynamics_cccv_complete(t, ...
%               I_chg, I_dch, I_current_cutoff, alpha, ...
%               curr.Ra, curr.Rb, curr.Qa, curr.Qb, za0, zb0, ...
%               f_ocv, Vmin, Vmax);

end

function curr = update_states_simple(prev, delta_time)
    % Simple form of the Q-R state update equation based on a constant
    % reaction rate constant k
    n = 0.5;

    k = 1;

    % Discrete state update equation for SEI growth
    curr.dQ = (k^(1/n) * delta_time + prev.dQ^(1/n))^n;
    
    % Ignore the resistance update for now; come back to this once we
    % figure out what we're doing for the capacity updates.
    curr.dR = prev.dR;

end

function curr = update_states_iss(prev, delta_time, Iss)
    % Q-R state update function based on steady-state current imbalance
    % values
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
    
    % Update rate constants for capacities and resistances
    gamma_q = 8;
    kq = gamma_q*abs(Iss);

    gamma_r = 0.00001;
    kr = gamma_r*abs(Iss);

    % Discrete state update equation for SEI growth
    curr.dQ = ( kq^(1/n) * delta_time + prev.dQ^(1/n) )^n;
    curr.dR = ( kr^(1/n) * delta_time + prev.dR^(1/n) )^n;
    
end