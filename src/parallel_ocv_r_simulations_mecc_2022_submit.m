function parallel_ocv_r_simulations_mecc_2022_submit()
    %% Andrew Weng
    %
    % Code used for MECC 2022 submission, April 22, 2022
    % Paper on parallel cell dynamics

    clear all; close all; clc
    set_default_plot_settings_manuscript()

    % Initialize model parameters
    Ra = 0.05; % ohms
    Qa = 5 * 3600; % As
    alpha = 1.2;
    za0 = 0.80;
    zb0 = 0.75;
    U0  = 3.0;
    current_target = -Qa / (3 * 3600);

    % Initialize simulation vectors
    t = linspace(0, 2.5*3600, 1e4)';
    I = current_target*ones(size(t)); % applied current (A) 

    %% Make the plots
    plot_dz_ss_sensitivity(Ra, current_target, alpha)
    plot_di_ss_sensitivity(Qa, Ra, current_target)
    plot_z_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)
    plot_i_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)
    
end

function plot_tau_sensitivity(Qa, Ra, alpha)

    %% Time constant plot
    rr = linspace(0.5, 1.5, 100);
    qq = linspace(0.5, 1.5, 101);

    for i = 1:numel(rr)
        for j = 1:numel(qq)

            Rb = Ra / rr(i);
            Qb = Qa / qq(j);

            tau(i, j) = (Ra+Rb)/alpha * (Qa*Qb/(Qa+Qb)) / 3600; % time constant in hours

        end
    end

    fh = figure('Position', [500 100 400*1.5 400*1.5]);
    [C,h] = contour(qq, rr, tau, 'LabelSpacing', 300, ...
                                 'ShowText', 'on', 'LineColor', 'k');
    clabel(C, h, 'FontSize', 24, 'Color', 'k', 'FontName', 'Times New Roman')
    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])
    yline(1, 'LineStyle', ':')
    xline(1, 'LineStyle', ':')

    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\tau$ (hrs)', 'Interpreter', 'Latex')
    saveas(fh, 'figs/fig_time_constant.png')

end

function plot_dz_ss_sensitivity(Ra, current_target, alpha)

    set_default_plot_settings_manuscript()

    rr = linspace(0.5, 1.5, 101);
    qq = linspace(0.5, 1.5, 100);

    % Note: this can be vectorized
    for i = 1:numel(rr)
        for j = 1:numel(qq)
            Rb = Ra / rr(i);
            dz_ss(i, j) = Rb/alpha * (rr(i)*qq(j) - 1) / (qq(j) + 1) * current_target;
        end
    end

    % Highlight the special solution for zero delta zz.
    rr_special = rr;
    qq_special = 1./rr_special;

    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    line(rr_special, qq_special, ...
        'LineWidth', 5, ...
        'Color',' k', ...
        'LineStyle', ':', ...
        'Parent', gca, ...
        'DisplayName', '$q=r^{-1}$'); hold all;

    [C,h] = contour(qq, rr, dz_ss, 'LabelSpacing', 300, ...
                                 'ShowText', 'on', 'LineColor', 'k', ...
                                 'Parent', gca, ...
                                 'HandleVisibility', 'off');

    clabel(C, h, 'FontSize', 24, 'Color', 'k', 'FontName', 'Times New Roman')

    q_vec = [1, 0.8, 1];
    r_vec = [1, 1.25, 1.5];
    marker_vec = {'o', 's', '^'};

    for i = 1:numel(q_vec)

        q = q_vec(i);
        r = r_vec(i);

        line(q, r, ...
        'Marker', marker_vec{i},  ...
        'Color', [0 0.5 0], ...
        'MarkerFaceColor', [0 0.5 0], ...
        'LineStyle', 'none', ...
        'MarkerSize', 12, ...
        'DisplayName', sprintf('$(q,r)=(%g,%g)$', q, r));

    end

    lh = legend('show'); set(lh, 'Interpreter', 'Latex', ...
                                 'Location', 'SouthWest')

    box on;

    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])
    yline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xlim([min(rr), max(rr)])
    ylim([min(qq), max(qq)])

    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\Delta z_{ss}$', 'Interpreter', 'Latex')

    saveas(fh, 'figs/fig_dz_ss_sensitivity.png')

end

function plot_di_ss_sensitivity(Qa, Ra, current_target)

    set_default_plot_settings_manuscript()

    %% Delta I_ss plot
    rr = linspace(0.5, 1.5, 101);
    qq = linspace(0.5, 1.5, 100);

    % Note: this can be vectorized
    for i = 1:numel(rr)
        for j = 1:numel(qq)
            Rb = Ra / rr(i);
            Qb = Qa / qq(j);
            di_ss(i, j) = current_target * ...
                (  2*(Ra*Qa - Rb*Qb)/((Ra+Rb)*(Qa+Qb)) - (Ra-Rb)/(Ra+Rb) );
        end
    end
   
    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    line(ones(size(qq)), qq, ...
        'LineWidth', 5, ...
        'Color',' k', ...
        'LineStyle', ':', ...
        'Parent', gca, ...
        'DisplayName', '$q=1$'); hold all;
    
    [C,h] = contour(qq, rr, di_ss, 'LabelSpacing', 300, ...
                                    'ShowText', 'on', 'LineColor', 'k', ...
                                    'Parent', gca, ...
                                    'HandleVisibility', 'off');

    box on;

    q_vec = [1, 0.8, 1];
    r_vec = [1, 1.25, 1.5];
    marker_vec = {'o', 's', '^'};

    for i = 1:numel(q_vec)

        q = q_vec(i);
        r = r_vec(i);

        line(q, r, ...
        'Marker', marker_vec{i},  ...
        'Color', [0 0.5 0], ...
        'MarkerFaceColor', [0 0.5 0], ...
        'LineStyle', 'none', ...
        'MarkerSize', 12, ...
        'DisplayName', sprintf('$(q,r)=(%g,%g)$', q, r));

    end


    lh = legend('show'); set(lh, 'Interpreter', 'Latex', ...
                                 'Location', 'SouthWest')

    clabel(C, h, 'FontSize', 24, 'Color', 'k', ...
        'FontName', 'Times New Roman', 'Interpreter', 'Latex')
    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])
    yline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xlim([min(rr), max(rr)])
    ylim([min(qq), max(qq)])
    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\Delta I_{ss}$ (A)', 'Interpreter', 'Latex')

    saveas(fh, 'figs/fig_di_ss_sensitivity.png')

end

function plot_z_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)

    linestyles = {'-', '--', ':'};
    marker_vec = {'o', 's', '^'};
    ms = 9;

    q_vec = [1, 0.8, 1];
    r_vec = [1, 1.25, 1.5];

    fh = figure('Position', [500 100 400*1.5 500*1.5]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    ax1 = nexttile(th, 1); box on;
    ax1 = gca(); set(ax1,'XTickLabel',[])
    ylabel(ax1, '$\Delta z$', 'Interpreter', 'Latex')
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    
    ax2 = nexttile(th, 2);
    ax2 = gca(); box on;
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$z$', 'Interpreter', 'Latex')    
    
    ax2_inset = axes('Position', [0.63 0.15 0.30 0.15], 'Units', 'Normalized'); 
    box on; 
    ax2_inset.FontSize = 14;
        
    for i = 1:numel(q_vec)

        q = q_vec(i);
        r = r_vec(i);
        
        Qb = Qa/q;
        Rb = Ra/r;

        label = sprintf('$(q,r)=(%g,%g)$', q, r);
    
        % Do the calculations
        [za, zb] = solve_z_dynamics(t, I, alpha, Ra, Rb, Qa, Qb, za0, zb0);
        dz = za - zb;
        [Ia, Ib] = solve_branch_currents(I, dz, alpha, Ra, Rb);
        Vta = U0 + alpha.*za - Ia.*Ra;

        time_at_vmax = interp1(Vta, t, 4.2);
    
        idx = find(t < time_at_vmax);

        this_t = t(idx)/3600;
        this_za = za(idx);
        this_zb = zb(idx);
        this_dz = this_za - this_zb;

        line(this_t, this_dz, 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Color', 'k', 'DisplayName', label, 'Parent', ax1)
        
        line(this_t, this_za, 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
        line(this_t, this_zb, 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
       
        line(t(idx)/3600, za(idx), 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2_inset)
        line(t(idx)/3600, zb(idx), 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2_inset)
  
        line(this_t(end), this_dz(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax1, 'HandleVisibility', 'off')
        line(this_t(end), this_za(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
        line(this_t(end), this_zb(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
        line(this_t(end), this_za(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax2_inset, 'HandleVisibility', 'off')
        line(this_t(end), this_zb(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax2_inset, 'HandleVisibility', 'off')
              
    end

    line(1000, 1000, ...
        'Marker', 'o', 'MarkerSize', ms, 'Color', [0 0.5 0], 'LineStyle', 'none', ...
        'MarkerFaceColor', [0 0.5 0], 'DisplayName', '$V_t=4.2V$', 'Parent', ax1)
    
    xlim(ax1, [0, 1.5])
    ylim(ax1, [-0.02 0.06])
    ylim(ax2, [0.7 1.02])
    linkaxes([ax1 ax2], 'x')
    xlim(ax2_inset, [1.05, 1.4])
    ylim(ax2_inset, [0.96, 0.98])

    legend(ax1, 'show', 'Location', 'NorthEast', 'Interpreter', 'Latex')
    legend(ax2, {'$z_a$', '$z_b$'}, 'Location', 'NorthWest', 'Interpreter', 'Latex')

    saveas(fh, sprintf('figs/fig_z_dynamics_timeseries.png'))

end

function plot_i_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)

    linestyles = {'-', '--', ':'};
    marker_vec = {'o', 's', '^'};
    ms = 9;

    q_vec = [1, 0.8, 1];
    r_vec = [1, 1.25, 1.5];

    fh = figure('Position', [500 100 400*1.5 500*1.5]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    ax1 = nexttile(th, 1); box on;
    ax1 = gca(); set(ax1,'XTickLabel',[])
    ylabel(ax1, '$\Delta I$ (A)', 'Interpreter', 'Latex')
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    
    ax2 = nexttile(th, 2);
    ax2 = gca(); box on;
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')    
        
    for i = 1:numel(q_vec)

        q = q_vec(i);
        r = r_vec(i);
        
        Qb = Qa/q;
        Rb = Ra/r;

        label = sprintf('$(q,r)=(%g,%g)$', q, r);
    
        % Do the calculations
        [za, zb] = solve_z_dynamics(t, I, alpha, Ra, Rb, Qa, Qb, za0, zb0);
        dz = za - zb;
        [Ia, Ib] = solve_branch_currents(I, dz, alpha, Ra, Rb);
        Vta = U0 + alpha.*za - Ia.*Ra;

        time_at_vmax = interp1(Vta, t, 4.2);
    
        idx = find(t < time_at_vmax);

        this_t = t(idx)/3600;
        this_ia = Ia(idx);
        this_ib = Ib(idx);
        this_di = this_ia - this_ib;

        line(this_t, this_di, 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Color', 'k', 'DisplayName', label, 'Parent', ax1)
        line(this_t, this_ia, 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
        line(this_t, this_ib, 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)

        line(this_t(end), this_di(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax1, 'HandleVisibility', 'off')
        line(this_t(end), this_ia(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
        line(this_t(end), this_ib(end), 'Marker', marker_vec{i}, 'MarkerSize', ms, 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
 
    end

    line(1000, 1000, ...
        'Marker', 'o', 'MarkerSize', ms, 'Color', [0 0.5 0], 'LineStyle', 'none', ...
        'MarkerFaceColor', [0 0.5 0], 'DisplayName', '$V_t=4.2V$', 'Parent', ax1)    

    xlim(ax1, [0, 1.5])
    ylim(ax1, [-0.1 2.0])
    ylim(ax2, [-1.7 0.2])
    linkaxes([ax1 ax2], 'x')

    legend(ax1, 'show', 'Location', 'NorthEast', 'Interpreter', 'Latex')
    legend(ax2, {'$I_a$', '$I_b$'}, 'Location', 'NorthEast', 'Interpreter', 'Latex')

    saveas(fh, sprintf('figs/fig_i_dynamics_timeseries.png'))

end