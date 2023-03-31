function fig_affine_dynamics()
% Affine dynamics illustration figures

    set_default_plot_settings();

    % Initialize model parameters
    Ra  = 0.15; % ohms
    Qa  = 3.00 * 3600; % As
    za0 = 0.20;
    zb0 = 0.30;
    q   = 0.7;
    r   = 1.1;

    Vmax = 4.2;
    Vmin = 3.0;
    U0 = 3.0;

    alpha = Vmax - Vmin;

    Qb = Qa/q;
    Rb = Ra/r;

    current_target = -Qa / (1 * 3600);

    %% Test the analytic solution
    ocv_lin = @(z) U0 + alpha * z;
    ocv_nonlin = load_ocv_fn('nmc-affine');

    % Initialize simulation parameters
    t = linspace(0, 10*3600, 1.0e6)';
    I_chg = +current_target*ones(size(t)); % applied current (A) 
    I_dch = -current_target*ones(size(t));
    I_cutoff = current_target/5;
    
    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);

    res_lsim2 = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, res_lsim.za(end), res_lsim.zb(end), ...
        ocv_lin, Vmin, Vmax);

    res_lsim3 = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, alpha, Ra, Rb, Qa, Qb, res_lsim2.za(end), res_lsim2.zb(end), ...
        ocv_lin, Vmin, Vmax);


    plot_results_default(res_lsim, ocv_lin)
    plot_results_imbalance(res_lsim)

    % Plot limit cycles in the z_a, z_b diagram
    fh = figure('Position', [500 100 500 500]); box on; grid off;

    line(linspace(0, 1, 100), linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_1=z_2$')
    line(res_lsim.zb(1), res_lsim.za(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
    line(res_lsim.zb, res_lsim.za, 'LineStyle', '-', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '1st Cycle') 
    line(res_lsim2.zb, res_lsim2.za, 'LineStyle', '--', 'LineWidth', 3, 'Color', 'r', 'DisplayName', '2nd Cycle') 
    line(res_lsim3.zb, res_lsim3.za, 'LineStyle', ':', 'LineWidth', 2, 'Color', 'b', 'DisplayName', '3rd Cycle') 

%     line(Qb.*linspace(0, 1, 100), Qa.*linspace(0, 1, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$z_1=z_2$')
%     line(Qb.*res_lsim.zb(1), Qa.*res_lsim.za(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
%     line(Qb.*res_lsim.zb, Qa.*res_lsim.za, 'LineStyle', '-', 'LineWidth', 3, 'Color', 'k', 'DisplayName', '1st Cycle') 
%     line(Qb.*res_lsim2.zb, Qa.*res_lsim2.za, 'LineStyle', '--', 'LineWidth', 3, 'Color', 'r', 'DisplayName', '2nd Cycle') 
%     line(Qb.*res_lsim3.zb, Qa.*res_lsim3.za, 'LineStyle', ':', 'LineWidth', 2, 'Color', 'b', 'DisplayName', '3rd Cycle') 


    legend show
    ZLIMS = [min([min(res_lsim.za) min(res_lsim.zb)]) - 0.03 1];
    xlim(ZLIMS)
    ylim(ZLIMS)
    xlabel('$z_1$', 'Interpreter', 'latex')
    ylabel('$z_2$', 'Interpreter', 'latex')



    % Plot limit cycles in the z_a, z_b diagram
%     fh = figure('Position', [500 100 500 500]); box on; grid off;
% 
%     line(linspace(-2, 2, 100), linspace(-2, 2, 100), 'LineStyle', '--', 'Color', [0.5, 0.5, 0.5], 'DisplayName', '$I_1=I_2$')
%     line(res_lsim.Ib(1), res_lsim.Ia(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
%     line(res_lsim.Ib, res_lsim.Ia, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 3, 'Color', 'k', 'DisplayName', '1st Cycle') 
%     line(res_lsim2.Ib, res_lsim2.Ia, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 3, 'Color', 'r', 'DisplayName', '2nd Cycle') 
%     line(res_lsim3.Ib, res_lsim3.Ia, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 2, 'Color', 'b', 'DisplayName', '3rd Cycle') 
% 
%     legend show
% %     ZLIMS = [min([min(res_lsim.za) min(res_lsim.zb)]) - 0.03 1];
% %     xlim(ZLIMS)
% %     ylim(ZLIMS)
%     xlabel('$I_1$', 'Interpreter', 'latex')
%     ylabel('$I_2$', 'Interpreter', 'latex')
% 
    % Plot limit cycles in the z_a, z_b diagram
    fh = figure('Position', [500 100 500 500]); box on; grid off;
    line(res_lsim.Ia(1) - res_lsim.Ib(1), res_lsim.za(1) - res_lsim.zb(1), 'Marker', 'x', 'MarkerSize', 18, 'LineStyle', 'none', 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'I.C.')
    line(res_lsim.Ia - res_lsim.Ib, res_lsim.za - res_lsim.zb, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 3, 'Color', 'k', 'DisplayName', '1st Cycle') 
    line(res_lsim2.Ia - res_lsim2.Ib, res_lsim2.za - res_lsim2.zb, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 3, 'Color', 'r', 'DisplayName', '2nd Cycle') 
    line(res_lsim3.Ia - res_lsim3.Ib, res_lsim3.za - res_lsim3.zb, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 3, 'LineWidth', 2, 'Color', 'b', 'DisplayName', '3rd Cycle') 
    xlabel('$\Delta I$', 'Interpreter', 'latex')
    ylabel('$\Delta z$', 'Interpreter', 'latex')


end

function plot_results_default(res, ocv)

    %% Plot the results
    fh = figure('Position', [500 100 600 800]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    FONT = 20;
    TLIM = [0 3.5];
    ILIM = [min([min(res.Ia), min(res.Ib)]) - 0.2 max([max(res.Ia) max(res.Ib)]) + 0.6];
    ZLIM = [min([min(res.za), min(res.zb)]) - 0.03 max([max(res.za) max(res.zb)]) + 0.06];
    VLIM = [min([min(ocv(res.za)), min(ocv(res.zb))]) - 0.1 max([max(ocv(res.za)) max(ocv(res.zb))]) + 0.1];

    % Current plot
    ax1 = nexttile(th, 1); set(ax1, 'XTickLabel', []); box on; %ylim([-1.5 1.2])
    ylabel(ax1, '$I$ (A)', 'Interpreter', 'Latex')
    
    line(res.t./3600, res.Ib, 'Color', 'b', 'DisplayName', '$I_1$', 'LineWidth', 2)
    line(res.t./3600, res.Ia, 'Color', 'r', 'DisplayName', '$I_2$', 'LineWidth', 2)

    if isfield(res, 'Ia_ohmic')
%         line(res.t./3600, res.Ib_ohmic, 'Color', 'b', 'DisplayName', '$I_{1,\mathrm{ohmic}}$', 'LineStyle', ':')
        line(res.t./3600, res.Ib_rebal, 'Color', 'b', 'DisplayName', '$I_{1,\mathrm{rebal}}$', 'LineStyle', '--')
%         line(res.t./3600, res.Ia_ohmic, 'Color', 'r', 'DisplayName', '$I_{2,\mathrm{ohmic}}$', 'LineStyle', ':', 'LineWidth', 2)
        line(res.t./3600, res.Ia_rebal, 'Color', 'r', 'DisplayName', '$I_{2,\mathrm{rebal}}$', 'LineStyle', '--', 'LineWidth', 2)

    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    ylim(ILIM)
    lh = legend('show', 'Location', 'east', 'FontSize', 20);

    % SOC plot
    ax2 = nexttile(th, 2); box on; set(ax2,'XTickLabel',[]); %ylim([0.78 1.05])
    ylabel('$z$', 'Interpreter', 'Latex')
    line(res.t./3600, res.zb, 'Color', 'b', 'DisplayName', '$z_1$', 'LineWidth', 2)
    line(res.t./3600, res.za, 'Color', 'r', 'DisplayName', '$z_2$', 'LineWidth', 2)
    if isfield(res, 'za_ohmic')
%         line(res.t./3600, res.za_ohmic + res.za(1), 'Color', 'r', 'DisplayName', '$z_{a,\mathrm{ext}}$', 'LineStyle', ':', 'LineWidth', 2)
%         line(res.t./3600, res.za_rebal, 'Color', 'r', 'DisplayName', '$z_{a,\mathrm{int}}$', 'LineStyle', '--', 'LineWidth', 2)
%         line(res.t./3600, res.zb_ohmic + res.zb(1), 'Color', 'b', 'DisplayName', '$z_{b,\mathrm{ext}}$', 'LineStyle', ':')
%         line(res.t./3600, res.zb_rebal, 'Color', 'b', 'DisplayName', '$z_{b,\mathrm{int}}$', 'LineStyle', '--')
    end

    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    ylim(ZLIM)
    lh = legend('show'); set(lh, 'FontSize', 20)

    % Voltage plot
    ax3 = nexttile(th, 3); box on; %ylim([3.9, 4.25])
    line(res.t./3600, res.Vt, 'Color', 'k', 'LineWidth', 2, 'DisplayName', '$V_t$')
    line(res.t./3600, ocv(res.zb), 'LineWidth', 2, 'Color', 'b', 'DisplayName', '$U_1$')
    line(res.t./3600, ocv(res.za), 'LineWidth', 2, 'Color', 'r', 'DisplayName', '$U_2$')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off')
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off')
    lh = legend('show', 'Location', 'northeast', 'FontSize', 20);
    xlabel('Time (hrs)')
    ylim(VLIM)
    ylabel('$V$ (V)', 'Interpreter', 'Latex')

    linkaxes([ax1, ax2, ax3], 'x')
    xlim(TLIM)
    xlabel('$t$ (hrs)', 'Interpreter', 'Latex')


end

function plot_results_imbalance(res)
    %% Make plots for the imbalance dynamics

    %% Plot the results
    fh = figure('Position', [500 100 600 600]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    FONT = 20;
    ILIM = [min(res.Ia - res.Ib) - 0.3 max(res.Ia - res.Ib) + 0.3];

   % Current plot
    ax1 = nexttile(th, 1); box on; ylim(ILIM)
    ylabel(ax1, '$\Delta I$ (A)', 'Interpreter', 'Latex')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    yline(0, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    
    line(res.t./3600, res.Ia - res.Ib, 'Color', 'k', 'LineWidth', 3, 'DisplayName', '$\Delta I$')
    if isfield(res, 'Ia_ohmic')
        line(res.t./3600, res.Ia_ohmic - res.Ib_ohmic, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', '$\Delta I_{\mathrm{ohmic}}$')
        line(res.t./3600, res.Ia_rebal - res.Ib_rebal, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '$\Delta I_{\mathrm{rebal}}$')
    end
    
    lh = legend('show', 'Location', 'best', 'FontSize', FONT);

    % SOC plot
    ax2 = nexttile(th, 2); box on; set(ax1,'XTickLabel',[]);
    ylabel('$\Delta z$', 'Interpreter', 'Latex')
    xline(res.t_chg_cc./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    xline(res.t_chg_cv./3600, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    yline(0, 'LineStyle', ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax2)
    
    line(res.t./3600, res.za - res.zb, 'Color', 'k', 'LineWidth', 3, 'DisplayName', '$\Delta z$')
    if isfield(res, 'za_ohmic')
%         line(res.t./3600, res.za_ohmic - res.zb_ohmic, 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', '$\Delta z_{\mathrm{ext}}$')
%         line(res.t./3600, res.za_rebal - res.zb_rebal, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', '$\Delta z_{\mathrm{int}}$')
    end

    lh = legend('Show', 'Location', 'Best', 'FontSize', FONT);
%     ylim([-0.02 0.02])
    xlabel('$t$ (hrs)', 'Interpreter', 'Latex')

    linkaxes([ax1, ax2], 'x')
    xlim([0 3.5])

end