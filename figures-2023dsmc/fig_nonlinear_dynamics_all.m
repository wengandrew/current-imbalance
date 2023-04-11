function fig_nonlinear_dynamics_all()
    % Run a full charge CCCV and discharge CC simulation

    set_default_plot_settings();

    % Initialize model parameters
    Ra  = 0.05; % ohms
    Qa  = 3.00 * 3600; % As
    za0 = 0.40;
    zb0 = 0.20;
    q   = 0.7;
    r   = 1.1;

    Qb = Qa/q;
    Rb = Ra/r;

    current_target = -Qa / (1 * 3600);

    %% Test the analytic solution
    ocv_lin = load_ocv_fn('nmc-affine');
    ocv_nmc = load_ocv_fn('nmc');
    ocv_lfp = load_ocv_fn('lfp');
    max_hours = 5;

    Vmax_affine = ocv_lin(1);
    Vmin_affine = ocv_lin(0);
    alpha = Vmax_affine - Vmin_affine;

    % Initialize simulation parameters
    t = linspace(0, max_hours*3600, 1.0e6)';

    I_chg = +current_target*ones(size(t)); % applied current (A)
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/20;

    res_lin = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin_affine, Vmax_affine);

    res_lin2 = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, res_lin.za(end), res_lin.zb(end), ocv_lin, Vmin_affine, Vmax_affine);

    res_nmc = run_discrete_time_simulation_multicycle(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nmc, 3.0, 4.2);

    res_nmc2 = run_discrete_time_simulation_complete(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, res_nmc.za(end), res_nmc.zb(end), ocv_nmc, 3.0, 4.2);

    res_lfp = run_discrete_time_simulation_multicycle(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_lfp, 3.0, 3.6);

    res_lfp2 = run_discrete_time_simulation_complete(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, res_lfp.za(end), res_lfp.zb(end), ocv_lfp, 3.0, 3.6);

    plot_results_imbalance(res_lin, res_nmc, res_lfp, max_hours)

    % ORBITS
    fh = figure('Position', [500 100 500 500]); box on; grid off;

    line(linspace(0, 1, 100), linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_a=z_b$')
    line(res_lin.zb(1), res_lin.za(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
    line(res_lin.zb, res_lin.za, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'k', 'DisplayName', 'Affine')
    line(res_nmc.zb, res_nmc.za, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'r', 'DisplayName', 'NMC/Gr')
    line(res_lfp.zb, res_lfp.za, 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b', 'DisplayName', 'LFP/Gr')
%     line(res_lin2.zb, res_lin2.za, 'LineStyle', '-', 'Marker', 'o', 'MarkerSize', 1, 'Color', 'k', 'DisplayName', 'Cycle 2 (Affine)')
%     line(res_nmc2.zb, res_nmc2.za, 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r', 'DisplayName', 'Cycle 2 (NMC/Gr)')
%     line(res_lfp2.zb, res_lfp2.za, 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b', 'DisplayName', 'Cycle 2 (LFP/Gr)')

    lh = legend('show'); set(lh, 'FontSize', 16)
    xlabel('$z_1$', 'Interpreter', 'Latex')
    ylabel('$z_2$', 'Interpreter', 'Latex')
    xlim([0 1])
    ylim([0 1])


end


function plot_results_imbalance(res_lin, res_nmc, res_lfp, max_hours)
    %% Make plots for the imbalance dynamics

    %% Plot the results
    fh = figure('Position', [500 100 600 600]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none');

    % Define time offsets so that the affine and the non-linear simulations
    % will share common points at the start and end of the CV hold phase.
    t0 = res_lin.t_chg_cc;
    t1 = res_lin.t_chg_cv;

    t2 = res_nmc.t(find(res_nmc.state==1, 1, 'first')); % end of charge cc
    t3 = res_nmc.t(find(res_nmc.state==2, 1, 'first')); % end of charge cv
    idxa = res_nmc.t < t3;
    idxb = res_nmc.t >= t3;

    t4 = res_lfp.t(find(res_lfp.state==1, 1, 'first')); % end of charge cc
    t5 = res_lfp.t(find(res_lfp.state==2, 1, 'first')); % end of charge cv
    idxc = res_lfp.t < t5;
    idxd = res_lfp.t >= t5;
    
    % Current plot
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);
    ylabel(ax1, '$\Delta I$ (A)', 'Interpreter', 'Latex')

    % Affine
    line(res_lin.t./3600 - (t0-t2)/3600, res_lin.Ia - res_lin.Ib, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Affine')

    % NMC
    line(res_nmc.t(idxa)./3600, res_nmc.Ia(idxa) - res_nmc.Ib(idxa), 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'NMC/Gr')
    line(res_nmc.t(idxb)./3600 + (t1-t3)/3600 - (t0-t2)/3600, res_nmc.Ia(idxb) - res_nmc.Ib(idxb), 'LineWidth', 2, 'Color', 'r', 'HandleVisibility', 'off')

    % LFP
    line(res_lfp.t(idxc)./3600 - (t4-t2)/3600, res_lfp.Ia(idxc) - res_lfp.Ib(idxc), 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'LFP/Gr')
    line(res_lfp.t(idxd)./3600 - (t0-t2)/3600 + (t1-t5)/3600, res_lfp.Ia(idxd) - res_lfp.Ib(idxd), 'LineWidth', 2, 'Color', 'b', 'HandleVisibility', 'off')
    
    % Decorators
    xline(res_lin.t_chg_cc./3600 - (t0-t2)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res_lin.t_chg_cv./3600 - (t0-t2)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('show', 'Location', 'best');
%     ylim([-0.4 +0.4])

    % SOC plot
    ax2 = nexttile(th, 2); box on; 
    xlabel(ax2, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel('$\Delta z$', 'Interpreter', 'Latex')

    line(res_lin.t./3600 - (t0-t2)/3600, res_lin.za - res_lin.zb, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Affine')
    line(res_nmc.t(idxa)./3600, res_nmc.za(idxa) - res_nmc.zb(idxa), 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'NMC/Gr')
    line(res_nmc.t(idxb)./3600 + (t1-t3)/3600 - (t0-t2)/3600, res_nmc.za(idxb) - res_nmc.zb(idxb), 'Color', 'r', 'LineWidth', 2, 'HandleVisibility', 'off')

    line(res_lfp.t(idxc)./3600 - (t4-t2)/3600, res_lfp.za(idxc) - res_lfp.zb(idxc), 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'LFP/Gr')
    line(res_lfp.t(idxd)./3600 - (t0-t2)/3600 + (t1-t5)/3600, res_lfp.za(idxd) - res_lfp.zb(idxd), 'LineWidth', 2, 'Color', 'b', 'HandleVisibility', 'off')
    
    xline(res_lin.t_chg_cc./3600 - (t0-t2)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res_lin.t_chg_cv./3600 - (t0-t2)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    lh = legend('Show', 'Location', 'Best');
%     ylim([-0.02 0.02])

    linkaxes([ax1, ax2], 'x')
%     xlim([0 max_hours])

end
