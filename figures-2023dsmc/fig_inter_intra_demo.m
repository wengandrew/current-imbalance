function fig_inter_intra_demo()
    % Demo the intra part of the inter-intra dynamics. Focus on fewer plots
    set_default_plot_settings();

    CHEMISTRY = 'nmc';

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

    % Override

%     Qa = 0.85206 * 3600;
%     Qb = 1.21724 * 3600;
%     Ra = 0.165465;
%     Rb = 0.158457;

    current_target = -3 * 3600 / (1 * 3600);

    %% Test the analytic solution
    ocv_lin = load_ocv_fn(affine_name);
    ocv_nonlin = load_ocv_fn(CHEMISTRY);
    max_hours = 20;

    % Initialize simulation parameters
    t = linspace(0, max_hours*3600, 1.0e6)';

    I_chg = +current_target*ones(size(t)); % applied current (A)
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/5;

    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);

    res_disc = run_discrete_time_simulation_multicycle(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_nonlin, Vmin, Vmax);

    plot_results_default(res_lsim, res_disc, ocv_lin, ocv_nonlin, max_hours)

end

function plot_results_default(resa, resb, ocv_lin, ocv_nonlin, max_hours)



    t0 = 0;
    t1 = 0;
    t2 = 0;
    t3 = 0;

    %% Plot the results
    fh = figure('Position', [500 100 400 600]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none');


    % Current plot
    ax1 = nexttile(th, 1); set(ax1, 'XTickLabel', []); box on;
    ylabel(ax1, '$I$ (A)', 'Interpreter', 'Latex')

    line(resb.t./3600, resb.Ib, 'Color', 'b', 'DisplayName', '$I_1$', 'LineWidth', 2)
    line(resb.t./3600, resb.Ia, 'Color', 'r', 'DisplayName', '$I_2$', 'LineWidth', 2)

    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    lh = legend('show', 'Location', 'east');

    % SOC plot
    ax2 = nexttile(th, 2); box on;
    ylabel('$z$', 'Interpreter', 'Latex')

    line(resb.t./3600, resb.zb, 'Color', 'b', 'DisplayName', '$z_1$', 'LineWidth', 2)
    line(resb.t./3600, resb.za, 'Color', 'r', 'DisplayName', '$z_2$', 'LineWidth', 2)

    xlabel('$t$ (hrs)', 'Interpreter', 'Latex')
    ylim([0.0 1.1])
    legend show

    linkaxes([ax1, ax2], 'x')
    xlim([0 max_hours])

end
