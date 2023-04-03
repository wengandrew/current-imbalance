function fig_nonlinear_dynamics_lfp()
    % Run a full charge CCCV and discharge CC simulation

    set_default_plot_settings();

    CHEMISTRY = 'lfp';

    % Initialize model parameters
    Ra  = 0.05; % ohms
    Qa  = 3.00 * 3600; % As
    za0 = 0.40;
    zb0 = 0.20;
    q   = 0.7;
    r   = 1.1;

    switch CHEMISTRY
        case 'lfp'
            Vmax = 3.6;
            U0 = 3.0;
            Vmin = 3.0;
            affine_name = 'lfp-affine';
        case 'nmc'
            Vmax = 4.2;
            Vmin = 3.0;
            U0 = 3.0;
            affine_name = 'nmc-affine';
        case 'nmc-umbl2022feb'
            Vmax = 4.2;
            Vmin = 3.0;
            U0 = 3.0;
            affine_name = 'nmc-affine';
    end

    Qb = Qa/q;
    Rb = Ra/r;

    keyboard

    current_target = -Qa / (1 * 3600);

    %% Test the analytic solution
    ocv_lin = load_ocv_fn(affine_name);
    ocv_nonlin = load_ocv_fn(CHEMISTRY);
    max_hours = 4.5;

    Vmax_affine = ocv_lin(1);
    Vmin_affine = ocv_lin(0);
    alpha = Vmax_affine - Vmin_affine;

    % Initialize simulation parameters
    t = linspace(0, max_hours*3600, 1.0e6)';

    I_chg = +current_target*ones(size(t)); % applied current (A)
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/20;

    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin_affine, Vmax_affine);

    res_disc = run_discrete_time_simulation_complete(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nonlin, Vmin, Vmax);

    plot_results_default(res_lsim, res_disc, ocv_lin, ocv_nonlin, max_hours)
    plot_results_imbalance(res_lsim, res_disc, max_hours)

end

function plot_results_default(resa, resb, ocv_lin, ocv_nonlin, max_hours)


    %% Plot the results
    fh = figure('Position', [500 100 600 800]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none');

    % Current plot
    ax1 = nexttile(th, 1); set(ax1, 'XTickLabel', []); box on;
    ylabel(ax1, '$I$ (A)', 'Interpreter', 'Latex')

    line(resb.t./3600, resb.Ib, 'Color', 'b', 'DisplayName', '$I_1$ (LFP/Gr)', 'LineWidth', 2)
    line(resb.t./3600, resb.Ia, 'Color', 'r', 'DisplayName', '$I_2$ (LFP/Gr)', 'LineWidth', 2)

    xline(resb.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(resb.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('show', 'Location', 'east');

    % SOC plot
    ax2 = nexttile(th, 2); box on; set(ax2,'XTickLabel',[]); %ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')

    line(resb.t./3600, resb.zb, 'Color', 'b', 'DisplayName', '$z_1$ (LFP/Gr)', 'LineWidth', 2)
    line(resb.t./3600, resb.za, 'Color', 'r', 'DisplayName', '$z_2$ (LFP/Gr)', 'LineWidth', 2)

    xline(resb.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(resb.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    ylim([0.0 1.1])
    legend show

    % Voltage plot
    ax3 = nexttile(th, 3); box on;
    xlabel(ax3, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax3, '$V$ (V)', 'Interpreter', 'Latex')

    line(resb.t./3600, resb.Vt, 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$V_t$')
    line(resb.t./3600, ocv_nonlin(resb.zb), 'LineWidth', 2, 'Color', 'b', 'DisplayName', '$U_1$ (LFP/Gr)')
    line(resb.t./3600, ocv_nonlin(resb.za), 'LineWidth', 2, 'Color', 'r', 'DisplayName', '$U_2$ (LFP/Gr)')
 
    line(resb.t./3600, resb.Vt, 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
    xline(resb.t_chg_cc./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    xline(resb.t_chg_cv./3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')

    legend show
    ylim([2.9 3.75])
    xlim([0 5])

    linkaxes([ax1, ax2, ax3], 'x')
    xlim([0 max_hours])

end

function plot_results_imbalance(resa, resb, max_hours)
    %% Make plots for the imbalance dynamics

    %% Plot the results
    fh = figure('Position', [500 100 600 600]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none');

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

    % ORBITS
    fh = figure('Position', [500 100 800 400]); box on; grid off;

    subplot(1, 2, 1); box on; grid off

    idx_chgcc = find(resb.state == 0);
    idx_chgcv = find(resb.state == 1);
    idx_dchcc = find(resb.state == 2);
    line(linspace(0, 1, 100), linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_a=z_b$')
    line(resa.zb(1), resa.za(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
    line(resa.zb, resa.za, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 1, 'Color', 'k', 'DisplayName', '1st Cycle (Affine)')
    line(resb.zb, resb.za, 'LineWidth', 2,    'Color',   'r', 'DisplayName', '1st Cycle (Non-linear)')
    legend show
    xlabel('$z_1$', 'Interpreter', 'Latex')
    ylabel('$z_2$', 'Interpreter', 'Latex')
    xlim([0 1])
    ylim([0 1])


    subplot(1, 2, 2); box on; grid off
    
    idx_chgcc = find(resb.state == 0);
    idx_chgcv = find(resb.state == 1);
    idx_dchcc = find(resb.state == 2);
    line(linspace(0, 1, 100), linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_a=z_b$')
    line(resa.zb(1), resa.za(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
    line(resa.zb, resa.za, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 1, 'Color', 'k', 'DisplayName', '1st Cycle (Affine)')
    line(resb.zb, resb.za, 'LineWidth', 2,    'Color',   'r', 'DisplayName', '1st Cycle (Non-linear)')
    xlabel('$z_1$', 'Interpreter', 'Latex')
    ylabel('$z_2$', 'Interpreter', 'Latex')
    xlim([0.85 1])
    ylim([0.85 1])

end
