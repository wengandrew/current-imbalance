function fig_affine_dynamics()
% Affine dynamics illustration figures

    set_default_plot_settings();

    % Initialize model parameters
    Ra  = 0.15; % ohms
    Qa  = 3.00 * 3600; % As
    za0 = 0.70;
    zb0 = 0.75;
    q   = 0.7;
    r   = 1.1;

    Vmax = 4.2;
    Vmin = 3.0;
    U0 = 3.0;

    alpha = Vmax - Vmin;

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

    res_lsim2 = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, res_lsim.za(end), res_lsim.zb(end), ...
        ocv_lin, Vmin, Vmax);


    plot_results_default(res_lsim)
    plot_results_imbalance(res_lsim)

    % Plot limit cycles in the z_a, z_b diagram
    fh = figure('Position', [500 100 400 400]); box on; grid off;
    line(linspace(0, 1, 100), linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_a=z_b$')
    line(res_lsim.za(1), res_lsim.zb(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Initial Condition')
    line(res_lsim.za, res_lsim.zb, 'LineStyle', '-', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '1st Cycle') 
    line(res_lsim2.za, res_lsim2.zb, 'LineStyle', ':', 'LineWidth', 2, 'Color', 'r', 'DisplayName', '2nd Cycle') 
    legend show
    xlim([0.6 1])
    ylim([0.6 1])
    xlabel('z_a')
    ylabel('z_b')

    % Plot the voltages
    fh = figure('Position', [500 100 400 400]); box on; grid off;
    line(res_lsim.t./3600, res_lsim.Vt, 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$V_t$')
    line(res_lsim.t./3600, ocv_lin(res_lsim.za), 'LineWidth', 2, 'Color', 'r', 'DisplayName', '$U_a$')
    line(res_lsim.t./3600, ocv_lin(res_lsim.zb), 'LineWidth', 2, 'Color', 'b', 'DisplayName', '$U_b$')
    xline(res_lsim.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    xline(res_lsim.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    xlabel('Time (hrs)')
    ylabel('Voltage (V)')

    legend show
    ylim([3.7 4.25])
    xlim([0 5])

    keyboard

end

function plot_results_default(res)

    %% Plot the results
    fh = figure('Position', [500 100 500 750]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'tight'); 
    
    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]); %ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.za, 'Color', 'r', 'DisplayName', '$z_a$', 'LineWidth', 2)
    line(res.t./3600, res.zb, 'Color', 'b', 'DisplayName', '$z_b$', 'LineWidth', 2)
    if isfield(res, 'za_ohmic')
%         line(res.t./3600, res.za_ohmic, 'Color', 'r', 'DisplayName', '$z_{a,\mathrm{ext}}$', 'LineStyle', ':', 'LineWidth', 2)
        line(res.t./3600, res.za_rebal, 'Color', 'r', 'DisplayName', '$z_{a,\mathrm{int}}$', 'LineStyle', '--', 'LineWidth', 2)
%         line(res.t./3600, res.zb_ohmic, 'Color', 'b', 'DisplayName', '$z_{b,\mathrm{ext}}$', 'LineStyle', ':')
        line(res.t./3600, res.zb_rebal, 'Color', 'b', 'DisplayName', '$z_{b,\mathrm{int}}$', 'LineStyle', '--')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    ylim([0.68 1.07])
    lh = legend('show'); set(lh, 'FontSize', 24)
    
    % Current plot
    ax2 = nexttile(th, 2); box on; %ylim([-1.5 1.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.Ia, 'Color', 'r', 'DisplayName', '$I_a$', 'LineWidth', 2)
    line(res.t./3600, res.Ib, 'Color', 'b', 'DisplayName', '$I_b$', 'LineWidth', 2)
    if isfield(res, 'Ia_ohmic')
%         line(res.t./3600, res.Ia_ohmic, 'Color', 'r', 'DisplayName', '$I_{a,\mathrm{ext}}$', 'LineStyle', ':', 'LineWidth', 2)
        line(res.t./3600, res.Ia_rebal, 'Color', 'r', 'DisplayName', '$I_{a,\mathrm{int}}$', 'LineStyle', '--', 'LineWidth', 2)
%         line(res.t./3600, res.Ib_ohmic, 'Color', 'b', 'DisplayName', '$I_{b,\mathrm{ext}}$', 'LineStyle', ':')
        line(res.t./3600, res.Ib_rebal, 'Color', 'b', 'DisplayName', '$I_{b,\mathrm{int}}$', 'LineStyle', '--')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'east', 'FontSize', 24);

    linkaxes([ax1, ax2], 'x')
    xlim([0 5])
    xlabel('Time (hrs)')


end

function plot_results_imbalance(res)
    %% Make plots for the imbalance dynamics

    %% Plot the results
    fh = figure('Position', [500 100 500 750]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'tight'); 
    
    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);
    ylabel('$\Delta z$', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.za - res.zb, 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\Delta z$')
    if isfield(res, 'za_ohmic')
        line(res.t./3600, res.za_ohmic - res.zb_ohmic, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', '$\Delta z_{\mathrm{ext}}$')
        line(res.t./3600, res.za_rebal - res.zb_rebal, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '$\Delta z_{\mathrm{int}}$')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('Show', 'Location', 'Best', 'FontSize', 24);
%     ylim([-0.02 0.02])

    % Current plot
    ax2 = nexttile(th, 2); box on; ylim([-0.4 0.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$\Delta I$ (A)', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.Ia - res.Ib, 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\Delta I$')
    if isfield(res, 'Ia_ohmic')
        line(res.t./3600, res.Ia_ohmic - res.Ib_ohmic, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', '$\Delta I_{\mathrm{ext}}$')
        line(res.t./3600, res.Ia_rebal - res.Ib_rebal, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '$\Delta I_{\mathrm{int}}$')
    end
    
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'best', 'FontSize', 24);
    xlabel('Time (hrs)')

    
    % Voltage plot
%     ax3 = nexttile(th, 3); box on; %ylim([3.9, 4.25])
%     line(res.t./3600, res.Vt, 'Color', 'k', 'DisplayName', '')
%     xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
%     xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
%     xlabel('Time (hrs)')
%     ylabel('Voltage (V)')
%     ylim([2.5 4.25])

    linkaxes([ax1, ax2], 'x')
    xlim([0 5])

end