function test_full_charge_discharge_simulations()
    % Verify the behavior of the code that runs the full charge-discharge
    % simulation

    set_default_plot_settings_manuscript()

    % Initialize model parameters
    Ra = 0.05; % ohms
    Qa = 5 * 3600; % As
    alpha = 1.2;
    za0 = 0.80;
    zb0 = 0.85;
    U0  = 3.0;
    Vmax = 4.2;
    Vmin = 3.0;

    q = 0.9;
    r = 1.5;
        
    Qb = Qa/q;
    Rb = Ra/r;

    current_target = -Qa / (3 * 3600);

    %% Test the analytic solution
    ocv_lin = @(z) U0 + alpha * z;

    % Initialize simulation vectors
    t = linspace(0, 10*3600, 1.0e5)';
    I_chg = +current_target*ones(size(t)); % applied current (A) 
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/20;
    out_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);


    %% Test the simulation-based solution
%     dt = 1;
%     out_disc = run_discrete_time_simulation_complete(dt, I_chg, I_dch, ...
%         I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_lin, Vmin, Vmax);

    %% Plot the results

    fh = figure('Position', [500 100 650 1000]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]); ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')
    
    line(out_lsim.t./3600, out_lsim.za, 'Color', 'r', 'DisplayName', '$z_a$')
    line(out_lsim.t./3600, out_lsim.zb, 'Color', 'b', 'DisplayName', '$z_b$', 'LineStyle', '-.')
%     line(out_disc.t./3600, out_disc.za, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
%     line(out_disc.t./3600, out_disc.zb, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')

    xline(out_lsim.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(out_lsim.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)

    legend show
    
    ax2 = nexttile(th, 2); box on; ylim([-1.5 1.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')
    
    line(out_lsim.t./3600, out_lsim.Ia, 'Color', 'r', 'DisplayName', '$I_a$')
    line(out_lsim.t./3600, out_lsim.Ib, 'Color', 'b', 'DisplayName', '$I_b$', 'LineStyle', '-.')
%     line(out_disc.t./3600, out_disc.Ia, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
%     line(out_disc.t./3600, out_disc.Ib, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
    
    xline(out_lsim.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(out_lsim.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'east');

%     ax3 = subplot(313);
%     line(out_lsim.t./3600, out_lsim.Vt, 'Color', 'k', 'DisplayName', 'Analytic')
%     line(out_disc.t./3600, out_disc.Vt, 'Color', 'k', 'Linestyle', ':', 'DisplayName', 'Simulated')
% 
%     xlabel('Time (hrs)')
%     ylabel('Voltage (V)')
%     legend show

    linkaxes([ax1, ax2], 'x')
    xlim([0 2.5])
%     xlim([0 3])

end