function fig_nonlin_highlights(chemistry, plot_type_int)
    % Run a full charge CCCV and discharge CC simulation
    %
    % Parameters
    % ----------
    % chemistry: 'lfp', 'nmc'
    % plot_type_int: (1) or (2)

    set_default_plot_settings();

    % Initialize model parameters
    Ra  = 0.15; % ohms
    Qa  = 3.00 * 3600; % As
    za0 = 0.40;
    zb0 = 0.20;
    q   = 0.7;

    switch plot_type_int
        case 1
            r = 1.1;
        case 2
            r = 1/q; % QR matching condition
    end

    switch chemistry
        case 'lfp'
            Vmax = 3.6;
            U0 = 3.0;
            Vmin = 3.0;
            label = '(LFP/Gr)';
            affine_name = 'lfp-affine';
        case 'nmc'
            Vmax = 4.2;
            Vmin = 3.31;
            U0 = 3.31;
            label = '(NMC/Gr)';
            affine_name = 'nmc-affine';
        case 'nmc-umbl2022feb'
            Vmax = 4.2;
            Vmin = 3.0;
            U0 = 3.0;
            label = '(NMC/Gr)';
            affine_name = 'nmc-affine';
    end

    Qb = Qa/q;
    Rb = Ra/r;

    p.Qa = Qa;
    p.Qb = Qb;
    p.Ra = Ra;
    p.Rb = Rb;

    kappa = 1 / (Vmax - Vmin) * (Ra*Qa - Rb*Qb) / (Qa + Qb);


    current_target = -Qa / (1 * 3600);

    %% Test the analytic solution
    ocv_lin = load_ocv_fn(affine_name);
    ocv_nonlin = load_ocv_fn(chemistry);
    max_hours = 15;

    % Initialize simulation parameters
    t = linspace(0, max_hours*3600, 1.0e6)';

    I_chg = +current_target*ones(size(t)); % applied current (A)
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/20;

    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);

    res_disc = run_discrete_time_simulation_multicycle(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nonlin, 3.0, Vmax);

    plot_results_default(res_lsim, res_disc, ocv_lin, ocv_nonlin, max_hours, p, za0-zb0, t, kappa, I_chg, label)
%     plot_results_imbalance(res_lsim, res_disc, max_hours)

end

function plot_results_default(resa, resb, ocv_lin, ocv_nonlin, max_hours, p, dz0, t, kappa, I, label)

    % Define time offsets so that the affine and the non-linear simulations
    % will share common points at the start and end of the CV hold phase.
    t0 = resa.t_chg_cc;
    t1 = resa.t_chg_cv;
    t3 = resb.t(find(resb.state==1, 1, 'first'));
    t2 = resb.t(find(resb.state==2, 1, 'first'));
    idxa = resb.t < t2;
    idxb = resb.t >= t2;

    %% Plot the results
    fh = figure('Position', [500 100 600 800]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none');


    % Current plot
    ax1 = nexttile(th, 1); set(ax1, 'XTickLabel', []); box on;
    ylabel(ax1, '$I$ (A)', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.Ib, 'Color', 'b', 'LineStyle', '--', 'DisplayName', '$I_1$ (affine)', 'LineWidth', 2)
    line(resa.t./3600 - (t0-t3)/3600, resa.Ia, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '$I_2$ (affine)', 'LineWidth', 2)

    line(resb.t(idxa)./3600, resb.Ib(idxa), 'Color', 'b', 'DisplayName', ['$I_1$ ' label], 'LineWidth', 2)
    line(resb.t(idxa)./3600, resb.Ia(idxa), 'Color', 'r', 'DisplayName', ['$I_2$ ' label]', 'LineWidth', 2)
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Ib(idxb), 'Color', 'b', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Ia(idxb), 'Color', 'r', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')


    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('show', 'Location', 'east');

    % SOC plot

    [condition, zbound_l2, zbound_linf, ibound_l2, ibound_linf] = ...
                solve_imbalance_bounds(ocv_nonlin, p, t, I, dz0);

    ax2 = nexttile(th, 2); box on; set(ax2,'XTickLabel',[]); %ylim([0.78 1.05])
    ylim([-0.23 0.23])
    ylabel('$\Delta z$', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.za - resa.zb, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta z$ (affine)')
    line(resb.t(idxa)./3600, resb.za(idxa) - resb.zb(idxa), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\Delta z$ (non-linear)')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.za(idxb) - resb.zb(idxb), 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off')

    yline(kappa*I(1), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Bound (affine, $\kappa$I)')
    line(t./3600,zbound_linf, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Bound (nonlinear)')
    yline(-kappa*I(1), 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility', 'off')
    line(t./3600,-zbound_linf, 'Color', 'r', 'LineWidth', 2, 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('Show', 'Location', 'Best');

    ylim([-0.5 0.5])
%     ylabel('$z$', 'Interpreter', 'Latex')
% 
%     line(resa.t./3600 - (t0-t3)/3600, resa.zb, 'Color', 'b', 'LineStyle', '--', 'DisplayName', '$z_1$ (affine)', 'LineWidth', 2)
%     line(resa.t./3600 - (t0-t3)/3600, resa.za, 'Color', 'r', 'LineStyle', '--', 'DisplayName', '$z_2$ (affine)', 'LineWidth', 2)
% 
% 
%     line(resb.t(idxa)./3600, resb.zb(idxa), 'Color', 'b', 'DisplayName', '$z_1$ (NMC/Gr)', 'LineWidth', 2)
%     line(resb.t(idxa)./3600, resb.za(idxa), 'Color', 'r', 'DisplayName', '$z_2$ (NMC/Gr)', 'LineWidth', 2)
% 
%     line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.za(idxb), 'Color', 'r', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')
%     line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.zb(idxb), 'Color', 'b', 'DisplayName', '', 'LineWidth', 2, 'HandleVisibility', 'off')
% 
%     xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
%     xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
%     ylim([0.0 1.1])
%     legend show

    % Voltage plot
    ax3 = nexttile(th, 3); box on;
    xlabel(ax3, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax3, '$V$ (V)', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.Vt, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$V_t$ (affine)')
    line(resa.t./3600 - (t0-t3)/3600, ocv_lin(resa.zb), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'b', 'DisplayName', '$U_1$ (affine)')
    line(resa.t./3600 - (t0-t3)/3600, ocv_lin(resa.za), 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r', 'DisplayName', '$U_2$ (affine)')
    line(resb.t(idxa)./3600, resb.Vt(idxa), 'Color', 'k', 'LineWidth', 2, 'DisplayName', ['$V_t$ ' label])
    line(resb.t(idxa)./3600, ocv_nonlin(resb.zb(idxa)), 'LineWidth', 2, 'Color', 'b', 'DisplayName', ['$U_1$ ' label])
    line(resb.t(idxa)./3600, ocv_nonlin(resb.za(idxa)), 'LineWidth', 2, 'Color', 'r', 'DisplayName', ['$U_2$ ' label])

    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Vt(idxb), 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, ocv_nonlin(resb.za(idxb)), 'LineWidth', 2, 'Color', 'r', 'HandleVisibility', 'off')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, ocv_nonlin(resb.zb(idxb)), 'LineWidth', 2, 'Color', 'b', 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')

    legend show
    ylim([2.9 4.25])
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

    % Current plot
    ax1 = nexttile(th, 1); box on; %ylim([-1.5 1.2])
    xlabel(ax1, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax1, '$\Delta I$ (A)', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.Ia - resa.Ib, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta I$ (affine)')
    line(resb.t(idxa)./3600, resb.Ia(idxa) - resb.Ib(idxa), 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta I$ (non-linear)')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.Ia(idxb) - resb.Ib(idxb), 'LineWidth', 2, 'Color', 'k', 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('show', 'Location', 'best');
%     ylim([-0.4 +0.4])


    % SOC plot
    ax2 = nexttile(th, 2); box on; set(ax2,'XTickLabel',[]);
    ylabel('$\Delta z$', 'Interpreter', 'Latex')

    line(resa.t./3600 - (t0-t3)/3600, resa.za - resa.zb, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', '$\Delta z$ (affine)')
    line(resb.t(idxa)./3600, resb.za(idxa) - resb.zb(idxa), 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$\Delta z$ (non-linear)')
    line(resb.t(idxb)./3600 + (t1-t2)/3600 - (t0-t3)/3600, resb.za(idxb) - resb.zb(idxb), 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off')

    xline(resa.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(resa.t_chg_cv./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('Show', 'Location', 'Best');
%     ylim([-0.02 0.02])

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
