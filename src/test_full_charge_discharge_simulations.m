function test_full_charge_discharge_simulations()
    % Verify the behavior of the code that runs the full charge-discharge
    % simulation

    set_default_plot_settings_manuscript()

    CHEMISTRY = 'nmc';

    % Initialize model parameters
    Ra = 0.05; % ohms
    Qa = 5 * 3600; % As
    za0 = 0.00;
    zb0 = 0.05;
    U0  = 3.0;

    switch CHEMISTRY
        case 'lfp'
            Vmax = 3.6;
        case 'nmc'
            Vmax = 4.2;
    end

    Vmin = 3.0;

    alpha = Vmax - Vmin;

    q = 0.9;
    r = 1.5;
        
    Qb = Qa/q;
    Rb = Ra/r;

    current_target = -Qa / (3 * 3600);

    %% Test the analytic solution
    ocv_lin = @(z) U0 + alpha * z;
    ocv_nonlin = load_ocv_fn(CHEMISTRY);

    % Initialize simulation parameters
    t = linspace(0, 10*3600, 1.0e6)';
    I_chg = +current_target*ones(size(t)); % applied current (A) 
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/20;
    
    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);

%     res_disc = run_discrete_time_simulation_complete(I_chg, I_dch, ...
%         I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nonlin, Vmin, Vmax);

    plot_results_default(res_lsim)
    plot_results_default(res_disc)
    plot_results_imbalance(res_lsim)

end

function plot_results_default(res)

    %% Plot the results
    fh = figure('Position', [500 100 650 900]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]); %ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.za, 'Color', 'r', 'DisplayName', '$z_a$')
    line(res.t./3600, res.zb, 'Color', 'b', 'DisplayName', '$z_b$', 'LineStyle', '-.')

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)

    legend show
    
    % Current plot
    ax2 = nexttile(th, 2); box on; set(ax2, 'XTickLabel',[]); %ylim([-1.5 1.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.Ia, 'Color', 'r', 'DisplayName', '$I_a$')
    line(res.t./3600, res.Ib, 'Color', 'b', 'DisplayName', '$I_b$', 'LineStyle', '-.')

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'east');

    % Voltage plot
    ax3 = nexttile(th, 3); box on; %ylim([3.9, 4.25])
    line(res.t./3600, gradient(res.Vt)./gradient(res.t./3600), 'Color', 'k', 'DisplayName', '')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xlabel('Time (hrs)')
    ylabel('Voltage (dV/dt)')

    linkaxes([ax1, ax2, ax3], 'x')
    %     xlim([0 2.5])
    
    % Here is some kind of experimental phase diagram
    figure(); plot(gradient(res.Vt)./gradient(res.t./3600), ...
              res.Ia - res.Ib, 'LineStyle', 'none', 'Marker', 'o', 'markerSize', 3, 'MarkerFaceColor', 'k', 'Color', 'k')
    xlabel('dV/dt (1/h)'); 
    ylabel('$\Delta I$ (A)', 'Interpreter', 'Latex')

end

function plot_results_imbalance(res)
    %% Make plots for the imbalance dynamics

    %% Plot the results
    fh = figure('Position', [500 100 650 900]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);
    ylabel('$\Delta z$', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.za - res.zb, 'Color', 'k', 'DisplayName', '$\Delta z$')
    line(res.t./3600, res.za_ohmic - res.zb_ohmic, 'Color', 'r', 'DisplayName', '$\Delta z_{\mathrm{ohmic}}$')
    line(res.t./3600, res.za_rebal - res.zb_rebal, 'Color', 'b', 'DisplayName', '$\Delta z_{\mathrm{rebal}}$')

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('Show', 'Location', 'Best');
    
    % Current plot
    ax2 = nexttile(th, 2); box on; set(ax2, 'XTickLabel',[]); %ylim([-1.5 1.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$\Delta I$ (A)', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.Ia - res.Ib, 'Color', 'k', 'DisplayName', '$\Delta I$')
    line(res.t./3600, res.Ia_ohmic - res.Ib_ohmic, 'Color', 'r', 'DisplayName', '$\Delta I_{\mathrm{ohmic}}$')
    line(res.t./3600, res.Ia_rebal - res.Ib_rebal, 'Color', 'b', 'DisplayName', '$\Delta I_{\mathrm{rebal}}$')

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'best');

    % Voltage plot
    ax3 = nexttile(th, 3); box on; %ylim([3.9, 4.25])
    line(res.t./3600, res.Vt, 'Color', 'k', 'DisplayName', '')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xlabel('Time (hrs)')
    ylabel('Voltage (V)')

    linkaxes([ax1, ax2, ax3], 'x')
%     xlim([0 2.5])

end