function fig_nonlinear_dynamics_lfp()
    % Run a full charge CCCV and discharge CC simulation

    set_default_plot_settings();

    CHEMISTRY = 'lfp';

    % Initialize model parameters
    Ra  = 0.15; % ohms
    Qa  = 3.00 * 3600; % As
    za0 = 0.0;
    zb0 = 0.0;
    q   = 0.7;
    r   = 1.1;

    switch CHEMISTRY
        case 'lfp'
            Vmax = 3.6;
            U0 = 2.7;
            Vmin = 2.7;
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

    max_hours = 17;

    % Initialize simulation parameters
    t = linspace(0, max_hours*3600, 1.0e6)';


    I_chg = +current_target*ones(size(t)); % applied current (A)
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/20;

    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);

    res_disc = run_discrete_time_simulation_complete(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nonlin, Vmin, Vmax);

    plot_results_default(res_lsim, res_disc, ocv_lin, max_hours)
    plot_results_imbalance(res_lsim, res_disc, max_hours)

end

function plot_results_default(resa, resb, ocv_lin, max_hours)

    % Define time offsets so that the affine and the non-linear simulations
    % will share common points at the start and end of the CV hold phase.
    t0 = resa.t_chg_cc;
    t1 = resa.t_chg_cv;
    t3 = resb.t(find(resb.state==1, 1, 'first'));
    t2 = resb.t(find(resb.state==2, 1, 'first'));
    idxa = resb.t < t2;
    idxb = resb.t >= t2;

    %% Plot the results
    fh = figure('Position', [500 100 700 700]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'tight');

    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]); %ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.za, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '$z_a$ (affine)', 'LineWidth', 2)
    line(resa.t./3600 - (t0-t3)/3600, resa.zb, 'Color', 'b', 'LineStyle', '--', 'DisplayName', '$z_b$ (affine)', 'LineWidth', 2)

    line(resb.t(idxa)./3600, resb.za(idxa), 'Color', 'r', 'DisplayName', '$z_a$ (non-linear)', 'LineWidth', 2)
    line(resb.t(idxa)./3600, resb.zb(idxa), 'Color', 'b', 'DisplayName', '$z_b$ (non-linear)', 'LineWidth', 2)
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.za(idxb), 'Color', 'r', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.zb(idxb), 'Color', 'b', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    ylim([0.0 1.09])
    legend show

    % Current plot
    ax2 = nexttile(th, 2); set(ax2, 'XTickLabel', []); box on;
    ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.Ia, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '$I_a$ (affine)', 'LineWidth', 2)
    line(resa.t./3600 - (t0-t3)/3600, resa.Ib, 'Color', 'b', 'LineStyle', '--', 'DisplayName', '$I_b$ (affine)', 'LineWidth', 2)
    line(resb.t(idxa)./3600, resb.Ia(idxa), 'Color', 'r', 'DisplayName', '$I_a$ (non-linear)', 'LineWidth', 2)
    line(resb.t(idxa)./3600, resb.Ib(idxa), 'Color', 'b', 'DisplayName', '$I_b$ (non-linear)', 'LineWidth', 2)
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Ia(idxb), 'Color', 'r', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Ib(idxb), 'Color', 'b', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'east');

    % Voltage plot
    ax3 = nexttile(th, 3); box on;
    xlabel(ax3, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax3, '$V$ (V)', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.Vt, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$V_t$ (affine)')
    line(resa.t./3600 - (t0-t3)/3600, ocv_lin(resa.za), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r', 'DisplayName', '$U_a$ (affine)')
    line(resa.t./3600 - (t0-t3)/3600, ocv_lin(resa.zb), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'b', 'DisplayName', '$U_b$ (affine)')
    line(resb.t(idxa)./3600, resb.Vt(idxa), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$V_t$ (non-linear)')
    line(resb.t(idxa)./3600, ocv_lin(resb.za(idxa)), 'LineWidth', 2, 'Color', 'r', 'DisplayName', '$U_a$ (non-linear)')
    line(resb.t(idxa)./3600, ocv_lin(resb.zb(idxa)), 'LineWidth', 2, 'Color', 'b', 'DisplayName', '$U_b$ (non-linear)')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Vt(idxb), 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, ocv_lin(resb.za(idxb)), 'LineWidth', 2, 'Color', 'r', 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, ocv_lin(resb.zb(idxb)), 'LineWidth', 2, 'Color', 'b', 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')

    legend show
    ylim([2.9 3.65])
    xlim([0 5])

    linkaxes([ax1, ax2, ax3], 'x')
    xlim([0 max_hours])

end

function plot_results_imbalance(resa, resb, max_hours)
    %% Make plots for the imbalance dynamics

    %% Plot the results
    fh = figure('Position', [500 100 700 700]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'tight');

    % Define time offsets so that the affine and the non-linear simulations
    % will share common points at the start and end of the CV hold phase.
    t0 = resa.t_chg_cc;
    t1 = resa.t_chg_cv;
    t3 = resb.t(find(resb.state==1, 1, 'first'));
    t2 = resb.t(find(resb.state==2, 1, 'first'));
    idxa = resb.t < t2;
    idxb = resb.t >= t2;

    % SOC plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);
    ylabel('$\Delta z$', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.za - resa.zb, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta z$ (affine)')
    line(resb.t(idxa)./3600, resb.za(idxa) - resb.zb(idxa), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\Delta z$ (non-linear)')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.za(idxb) - resb.zb(idxb), 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)

    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('Show', 'Location', 'Best');
%     ylim([-0.02 0.02])

    % Current plot
    ax2 = nexttile(th, 2); box on; %ylim([-1.5 1.2])
    xlabel(ax2, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$\Delta I$ (A)', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.Ia - resa.Ib, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta I$ (affine)')
    line(resb.t(idxa)./3600, resb.Ia(idxa) - resb.Ib(idxa), 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta I$ (non-linear)')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Ia(idxb) - resb.Ib(idxb), 'LineWidth', 2, 'Color', 'k', 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('show', 'Location', 'best');
%     ylim([-0.4 +0.4])

    linkaxes([ax1, ax2], 'x')
    xlim([0 max_hours])

    % SOC phase portrait
    ax3 = nexttile(th, 3); box on;
    idx_chgcc = find(resb.state == 0);
    idx_chgcv = find(resb.state == 1);
    idx_dchcc = find(resb.state == 2);
    line(linspace(0, 1, 100), linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_a=z_b$')
    line(resa.za(1), resa.zb(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Initial Condition')
    line(resa.za, resa.zb, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 1, 'Color', 'k', 'DisplayName', 'Affine')
    line(resb.za, resb.zb, 'LineWidth', 2,    'Color',   'r', 'DisplayName', 'Non-linear')
    legend show
    xlabel('z_a')
    ylabel('z_b')
    xlim([0 1])
    ylim([0 1])

end
