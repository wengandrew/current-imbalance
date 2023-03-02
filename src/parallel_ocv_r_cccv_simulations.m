function parallel_ocv_r_cccv_simulations()
    % Verify the behavior of the code that runs CC-CV simulations

    set_default_plot_settings_manuscript()

    analysis_type = 'lfp';

    switch analysis_type
        case 'lfp'
            alpha = 0.6;
            Vmax = 3.6;
            ocv = load_ocv_fn('lfp');
        case 'nmc'
            alpha = 1.2;
            Vmax = 4.2;
            ocv = load_ocv_fn('nmc');
    end

    % Initialize model parameters
    Ra = 0.05; % ohms
    Qa = 5 * 3600; % As
    alpha_nmc = 1.2;
    alpha_lfp = 0.6;
    za0 = 0.85;
    zb0 = 0.75;
    U0  = 3.0;
    Vmin = 3.0;

    q = 0.6;
    r = 1.5;
        
    Qb = Qa/q;
    Rb = Ra/r;
    
    current_target = -Qa / (3 * 3600);

    
    %% Test the analytic solution
    % Initialize simulation vectors
    t = linspace(0, 6*3600, 5e4)';
    I = current_target*ones(size(t)); % applied current (A) 

    [tfinal, za, zb, Ia, Ib] = solve_z_dynamics_cccv(t, I, ...
        alpha, Ra, Rb, Qa, Qb, za0, zb0, Vmax, U0);


    %% Test the simulation based solution
    out = run_discrete_time_simulation_cccv(t, I, Qa, Qb, Ra, Rb, ...
            za0, zb0, ocv, Vmax);


    % See the results
    figure()
    ax1 = subplot(311);
    line(tfinal./3600, za, 'Color', 'b', 'DisplayName', 'Analytic')
    line(tfinal./3600, zb, 'Color', 'r', 'DisplayName', 'Analytic')
    line(out.t./3600, out.za, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
    line(out.t./3600, out.zb, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
    xlabel('Time (hrs)')
    ylabel('SOC')
    legend show

    ax2 = subplot(312);
    line(tfinal./3600, Ia, 'Color', 'b', 'DisplayName', 'Analytic')
    line(tfinal./3600, Ib, 'Color', 'r', 'DisplayName', 'Analytic')
    line(out.t./3600, out.Ia, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
    line(out.t./3600, out.Ib, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
    xlabel('Time (hrs)')
    ylabel('Current (A)')
    legend show

    ax3 = subplot(313);
    line(tfinal/3600, U0 + alpha*za - Ia*Ra, 'Color', 'k', 'DisplayName', 'Analytic')
    line(out.t./3600, out.Vt, 'Color', 'k', 'LineStyle', ':', 'DisplayName', 'Simulated')
    xlabel('Time (hrs)')
    ylabel('Voltage (V)')
    legend show

    linkaxes([ax1, ax2, ax3], 'x')
    xlim([0 3])

end