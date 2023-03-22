function parallel_ocv_r_full_charge_discharge_simulation()
    % Run a full charge CCCV and discharge CC simulation

    set_default_plot_settings()

    CHEMISTRY = 'nmc';

    % Initialize model parameters
    Ra  = 0.164; % ohms
    Qa  = 1.98 * 3600; % As
    za0 = 0.50;
    zb0 = 0.35;
    q   = 0.7;
    r   = 1.1;

    switch CHEMISTRY
        case 'lfp'
            Vmax = 3.6;
            U0 = 3.0;
            Vmin = 3.0;
        case 'nmc'
            Vmax = 4.2;
            Vmin = 3.0;
            U0 = 3.0;
        case 'nmc-umbl2022feb'
            Vmax = 4.2;
            Vmin = 3.0;
            U0 = 3.0;
    end

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

    res_disc = run_discrete_time_simulation_complete(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nonlin, Vmin, Vmax);

    plot_results_default(res_lsim)
    plot_results_default(res_disc)
    plot_results_imbalance(res_lsim)
    plot_results_imbalance(res_disc)
    plot_results_phase(res_disc)

    idx_chgcc = find(res_disc.state == 0);
    idx_chgcv = find(res_disc.state == 1);
    idx_dchcc = find(res_disc.state == 2);

    fh = figure('Position', [500 100 600 600]); box on; grid on;
    line(res_lsim.za, res_lsim.zb, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 1, 'Color', 'k', 'DisplayName', 'Analytic')
    line(res_disc.za(idx_chgcc), res_disc.zb(idx_chgcc), 'LineWidth', 2,    'Color',   'r', 'DisplayName', 'Nonlin CCC')
    line(res_disc.za(idx_chgcv), res_disc.zb(idx_chgcv), 'LineStyle', ':',  'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Nonlin CCV')
    line(res_disc.za(idx_dchcc), res_disc.zb(idx_dchcc), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Nonlin DCC')
    legend show
    xlabel('z_a')
    ylabel('z_b')

    fh = figure('Position', [500 100 600 600]); box on; grid on;
    line(res_lsim.Ia, res_lsim.Ib, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 1, 'Color', 'k', 'DisplayName', 'Analytic')
    line(res_disc.Ia(idx_chgcc), res_disc.Ib(idx_chgcc), 'LineWidth', 2,    'Color',   'r', 'DisplayName', 'Nonlin CCC')
    line(res_disc.Ia(idx_chgcv), res_disc.Ib(idx_chgcv), 'LineStyle', ':',  'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Nonlin CCV')
    line(res_disc.Ia(idx_dchcc), res_disc.Ib(idx_dchcc), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Nonlin DCC')
    legend show
    xlabel('I_a')
    ylabel('I_b')

    keyboard

end

function plot_results_phase(res)

    % Here is some kind of experimental phase diagram
    figure();

    plot(gradient(res.Vt)./gradient(res.t./3600), ...
              res.Ia - res.Ib, 'LineStyle', 'none', ...
              'Marker', 'o', 'markerSize', 3, ...
              'MarkerFaceColor', 'k', 'Color', 'k')

    xlabel('dV/dt (V/h)');
    ylabel('$\Delta I$ (A)', 'Interpreter', 'Latex')
    ylim([-0.25 0.25]);
    xlim([-1.5 1.5])

end

function plot_results_default(res)

    %% Plot the results
    fh = figure('Position', [500 100 650 900]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none');

    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]); %ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')

    line(res.t./3600, res.za, 'Color', 'r', 'DisplayName', '$z_a$', 'LineWidth', 2)
    line(res.t./3600, res.zb, 'Color', 'b', 'DisplayName', '$z_b$', 'LineWidth', 2)
    if isfield(res, 'za_ohmic')
        line(res.t./3600, res.za_ohmic, 'Color', 'r', 'DisplayName', '$z_{a,\mathrm{ohmic}}$', 'LineStyle', ':', 'LineWidth', 2)
        line(res.t./3600, res.za_rebal, 'Color', 'r', 'DisplayName', '$z_{a,\mathrm{rebal}}$', 'LineStyle', '--', 'LineWidth', 2)
        line(res.t./3600, res.zb_ohmic, 'Color', 'b', 'DisplayName', '$z_{b,\mathrm{ohmic}}$', 'LineStyle', ':')
        line(res.t./3600, res.zb_rebal, 'Color', 'b', 'DisplayName', '$z_{b,\mathrm{rebal}}$', 'LineStyle', '--')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    ylim([-1.01 1.01])
    legend show

    % Current plot
    ax2 = nexttile(th, 2); box on; set(ax2, 'XTickLabel',[]); %ylim([-1.5 1.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')

    line(res.t./3600, res.Ia, 'Color', 'r', 'DisplayName', '$I_a$', 'LineWidth', 2)
    line(res.t./3600, res.Ib, 'Color', 'b', 'DisplayName', '$I_b$', 'LineWidth', 2)
    if isfield(res, 'Ia_ohmic')
        line(res.t./3600, res.Ia_ohmic, 'Color', 'r', 'DisplayName', '$I_{a,\mathrm{ohmic}}$', 'LineStyle', ':', 'LineWidth', 2)
        line(res.t./3600, res.Ia_rebal, 'Color', 'r', 'DisplayName', '$I_{a,\mathrm{rebal}}$', 'LineStyle', '--', 'LineWidth', 2)
        line(res.t./3600, res.Ib_ohmic, 'Color', 'b', 'DisplayName', '$I_{b,\mathrm{ohmic}}$', 'LineStyle', ':')
        line(res.t./3600, res.Ib_rebal, 'Color', 'b', 'DisplayName', '$I_{b,\mathrm{rebal}}$', 'LineStyle', '--')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'east');

    % Voltage plot
    ax3 = nexttile(th, 3); box on; %ylim([3.9, 4.25])
    line(res.t./3600, res.Vt, 'Color', 'k', 'DisplayName', '')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xlabel('Time (hrs)')
    ylabel('Voltage (V)')
    legend(sprintf('$t_{cc}$ = %.2f hrs \n$t_{cv}$ = %.2f hrs', ...
           res.t_chg_cc./3600, res.t_chg_cv./3600), ...
        'Interpreter', 'latex')
    ylim([2.5 4.25])

    linkaxes([ax1, ax2, ax3], 'x')
    %     xlim([0 2.5])

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
    if isfield(res, 'za_ohmic')
        line(res.t./3600, res.za_ohmic - res.zb_ohmic, 'Color', 'r', 'DisplayName', '$\Delta z_{\mathrm{ohmic}}$')
        line(res.t./3600, res.za_rebal - res.zb_rebal, 'Color', 'b', 'DisplayName', '$\Delta z_{\mathrm{rebal}}$')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('Show', 'Location', 'Best');
    ylim([-0.02 0.02])

    % Current plot
    ax2 = nexttile(th, 2); box on; set(ax2, 'XTickLabel',[]); %ylim([-1.5 1.2])
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$\Delta I$ (A)', 'Interpreter', 'Latex')

    line(res.t./3600, res.Ia - res.Ib, 'Color', 'k', 'DisplayName', '$\Delta I$')
    if isfield(res, 'Ia_ohmic')
        line(res.t./3600, res.Ia_ohmic - res.Ib_ohmic, 'Color', 'r', 'DisplayName', '$\Delta I_{\mathrm{ohmic}}$')
        line(res.t./3600, res.Ia_rebal - res.Ib_rebal, 'Color', 'b', 'DisplayName', '$\Delta I_{\mathrm{rebal}}$')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'best');
    ylim([-0.4 +0.4])

    % Voltage plot
    ax3 = nexttile(th, 3); box on; %ylim([3.9, 4.25])
    line(res.t./3600, res.Vt, 'Color', 'k', 'DisplayName', '')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax3)
    xlabel('Time (hrs)')
    ylabel('Voltage (V)')
    ylim([2.5 4.25])

    linkaxes([ax1, ax2, ax3], 'x')
%     xlim([0 2.5])

end
