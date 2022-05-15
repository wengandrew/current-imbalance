function parallel_ocv_r_simulations()
    %% Andrew Weng
    %% 4/1/2022

    clear all; close all; clc
    set_default_plot_settings_manuscript()

    % Initialize model parameters
    Ra = 0.05; % ohms
    Qa = 5 * 3600; % As
    alpha = 1.2;
    za0 = 0.8;
    zb0 = 0.75;
    U0 = 3.0;
    current_target = - Qa / (3 * 3600);

    % Initialize simulation vectors
    t = linspace(0, 2.5*3600, 1e4)';
    I = current_target*ones(size(t)); % applied current (A) 

    %% Make the plots

    plot_tau_sensitivity(Qa, Ra, alpha)

    % Plot dz and dI sensitivity in one figure

    plot_dz_ss_sensitivity(Ra, current_target, alpha)
    plot_di_ss_sensitivity(Qa, Ra, current_target)
    plot_di_ss_sensitivity_2(Qa)

%     plot_di_sensitivity(Qa, Ra, alpha)

    plot_z_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)
    plot_i_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)
    
    %% Tau contributions
    plot_tau_contributions()
end

function plot_tau_contributions()

    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    box on;
    x = linspace(0.5, 1.5, 100);
    rr = @(x) (x + 1) ./ x; 
    qq = @(x) 1 ./ (x + 1);

    line(x, rr(x), 'Color', 'k', 'LineWidth', 2)
    line(x, qq(x), 'Color', 'k', 'LineWidth', 2, 'LineStyle', ':')

    xlabel('$r$ or $q$', 'Interpreter', 'Latex')
    legend('$(r+1)/r$', '$1/(q+1)$', 'Interpreter', 'Latex')
   
end

function [Ia, Ib] = branch_currents(I, dz, alpha, Ra, Rb)

    Ia = (+alpha*dz + I*Rb) / (Ra + Rb);
    Ib = (-alpha*dz + I*Ra) / (Ra + Rb);

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
    saveas(fh, 'figures/fig_time_constant.png')

end

function plot_dz_ss_sensitivity(Ra, current_target, alpha)

    set_default_plot_settings_manuscript()

    rr = linspace(0.5, 1.5, 101);
    qq = linspace(0.5, 1.5, 100);

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
        'Color',' r', ...
        'LineStyle', ':', ...
        'Parent', gca); hold all;

    [C,h] = contour(qq, rr, dz_ss, 'LabelSpacing', 300, ...
                                 'ShowText', 'on', 'LineColor', 'k', ...
                                 'Parent', gca);

    clabel(C, h, 'FontSize', 24, 'Color', 'k', 'FontName', 'Times New Roman')


    box on;

    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])
    yline(1, 'LineStyle', ':')
    xline(1, 'LineStyle', ':')
    xlim([min(rr), max(rr)])
    ylim([min(qq), max(qq)])

    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\Delta z_{ss}$', 'Interpreter', 'Latex')

    saveas(fh, 'figures/fig_dz_ss_sensitivity.png')

end

function plot_di_ss_sensitivity(Qa, Ra, current_target)

    set_default_plot_settings_manuscript()

    %% Delta I_ss plot
    rr = linspace(0.5, 1.5, 101);
    qq = linspace(0.5, 1.5, 100);

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
        'Color',' r', ...
        'LineStyle', ':', ...
        'Parent', gca); hold all;
    
    [C,h] = contour(qq, rr, di_ss, 'LabelSpacing', 300, ...
                                    'ShowText', 'on', 'LineColor', 'k', ...
                                    'Parent', gca);

    box on;
    clabel(C, h, 'FontSize', 24, 'Color', 'k', 'FontName', 'Times New Roman', 'Interpreter', 'Latex')
    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])
    yline(1, 'LineStyle', ':')
    xline(1, 'LineStyle', ':')
    xlim([min(rr), max(rr)])
    ylim([min(qq), max(qq)])
    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\Delta I_{ss}$ (A)', 'Interpreter', 'Latex')

    saveas(fh, 'figures/fig_di_ss_sensitivity.png')

end

function plot_di_ss_sensitivity_2(Qa)

    set_default_plot_settings_manuscript()

    c_rate_vec = [-2, -1, 0, 1, 2];

    %% Delta I_ss plot
    qq = linspace(0.5, 1.5, 100);

    cmap = cool(numel(c_rate_vec));

    fh = figure('Position', [500 100 400*1.5 400*1.5]);
    box on;

    for i = 1:numel(c_rate_vec)
        
        current = c_rate_vec(i) * Qa / 3600;
        di_ss = (qq - 1) ./ (qq + 1) .* current;
        line(qq, di_ss, 'Color', cmap(i, :), ...
                        'DisplayName', sprintf('I = %.1fC', c_rate_vec(i)));
            
    end
    
    xlabel('$q$ ($Q_a/Q_b$)', 'Interpreter', 'Latex')
    ylabel('$\Delta I_{ss}$ (A)', 'Interpreter', 'Latex')
    lh = legend('show'); set(lh, 'Interpreter', 'Latex', 'Location', 'south', 'FontSize', 12)

    saveas(fh, 'figures/fig_di_ss_sensitivity_2.png')

end

function plot_di_sensitivity(Qa, Ra, alpha)

    set_default_plot_settings_manuscript()

    Rb = 1.5*Ra;
    dz_vec = -0.3:0.05:0.3;
    I_vec = linspace(-4*Qa/3600, 4*Qa/3600, 10);

    for i = 1:numel(dz_vec)
        for j = 1:numel(I_vec)
            
            dzz = dz_vec(i);
            ii = I_vec(j);
            dI(i, j) = dI_dynamics_zz(dzz, ii, Ra, Rb, alpha);

        end
    end

    fh = figure('Position', [500 100 400*1.5 400*1.5]);
    [C,h] = contour(I_vec, dz_vec, dI, 'LabelSpacing', 300, ...
                                        'ShowText', 'on', ...
                                        'LineColor', 'k');

    clabel(C, h, 'FontSize', 24, 'Color', 'k', 'FontName', 'Times New Roman')
    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])
    yline(0, 'LineStyle', ':')
    xline(0, 'LineStyle', ':')
    title(sprintf('$\\Delta I$ [A] \n ($R_a$ = %g m$\\Omega$, $R_b$ = %g m$\\Omega$)', ...
        Ra*1e3, Rb*1e3), 'Interpreter', 'Latex')
    ylabel('$\Delta z$', 'Interpreter', 'Latex')
    xlabel('$I$ (A)', 'Interpreter', 'Latex')
    saveas(fh, 'figures/fig_dI_dynamics.png')

end


function plot_z_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)

    linestyles = {'-', '--', ':'};

    q_vec = [1, 0.8, 1.0];    
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
        line(this_t(end), this_dz(end), 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineStyle', 'none', 'Parent', ax1, 'HandleVisibility', 'off')
        
        line(this_t, this_za, 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
        line(this_t, this_zb, 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
        line(this_t(end), this_za(end), 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
        line(this_t(end), this_zb(end), 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')

        line(t(idx)/3600, za(idx), 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2_inset)
        line(t(idx)/3600, zb(idx), 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2_inset)
        line(this_t(end), this_za(end), 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'LineStyle', 'none', 'Parent', ax2_inset, 'HandleVisibility', 'off')
        line(this_t(end), this_zb(end), 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'LineStyle', 'none', 'Parent', ax2_inset, 'HandleVisibility', 'off')
              
    end

    line(1000, 1000, 'Marker', 'o', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'DisplayName', '$V_t=4.2V$', 'Parent', ax1)
    
    xlim(ax1, [0, 1.5])
    ylim(ax1, [-0.03 0.06])
    ylim(ax2, [0.7 1.02])
    linkaxes([ax1 ax2], 'x')
    xlim(ax2_inset, [0.9, 1.4])
    ylim(ax2_inset, [0.95, 0.98])

    legend(ax1, 'show', 'Location', 'NorthEast', 'Interpreter', 'Latex')
    legend(ax2, {'$z_a$', '$z_b$'}, 'Location', 'NorthWest', 'Interpreter', 'Latex')

    saveas(fh, sprintf('figures/fig_z_dynamics_timeseries.png'))

end

function plot_i_timeseries(t, I, alpha, Ra, Qa, za0, zb0, U0)

    linestyles = {'-', '--', ':'};
    
    q_vec = [1, 0.8, 1.0];    
    r_vec = [1, 1.25, 1.5];

    fh = figure('Position', [500 100 400*1.5 500*1.5]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); 
    
    ax1 = nexttile(th, 1); box on;
    ax1 = gca(); set(ax1,'XTickLabel',[])
    ylabel(ax1, '$\Delta I$', 'Interpreter', 'Latex')
    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off', 'Parent', ax1)
    
    ax2 = nexttile(th, 2);
    ax2 = gca(); box on;
    xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$I$', 'Interpreter', 'Latex')    
        
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
        line(this_t(end), this_di(end), 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'k', 'LineStyle', 'none', 'Parent', ax1, 'HandleVisibility', 'off')
        
        line(this_t, this_ia, 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
        line(this_t, this_ib, 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
        line(this_t(end), this_ia(end), 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
        line(this_t(end), this_ib(end), 'Marker', 'o', 'Color', 'b', 'MarkerFaceColor', 'b', 'LineStyle', 'none', 'Parent', ax2, 'HandleVisibility', 'off')
              
    end

    line(1000, 1000, 'Marker', 'o', 'Color', 'k', 'LineStyle', 'none', 'MarkerFaceColor', 'k', 'DisplayName', '$V_t=4.2V$', 'Parent', ax1)
    
    xlim(ax1, [0, 1.5])
    ylim(ax1, [-0.3 1.6])
    ylim(ax2, [-1.7 0.2])
    linkaxes([ax1 ax2], 'x')

    legend(ax1, 'show', 'Location', 'NorthEast', 'Interpreter', 'Latex')
    legend(ax2, {'$I_a$', '$I_b$'}, 'Location', 'NorthEast', 'Interpreter', 'Latex')

    saveas(fh, sprintf('figures/fig_i_dynamics_timeseries.png'))

end

function plot_i_timeseries_old(t, I, alpha, Ra, Qa, za0, zb0, U0)

    %% Current Dynamics timeseries plots

    r_vec = [1, 1.5];
    linestyles = {'-', '--'};
    qb_vec = [Qa, 1.5*Qa];

    for iqb = 1:numel(qb_vec)

        Qb = qb_vec(iqb);
        capacity_string = sprintf('$Q_a$ = %g Ah, $Q_b$ = %g Ah', Qa/3600, Qb/3600);
    
        fh = figure('Position', [500 100 400*1.5 400*1.5]);
        th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'none'); 

        % Delta Current Dynamics
        nexttile(th, 1)
        
        ax1 = gca(); set(ax1,'XTickLabel',[])
        ylim(ax1, [-1, 6])

        title(capacity_string, 'Interpreter', 'Latex')
    
        for i = 1:numel(r_vec)

            Rb = r_vec(i)*Ra;
            dz = solve_dz_dynamics(t, I, alpha, Ra, Rb, Qa, Qb, za0 - zb0);
            [Ia, Ib] = solve_branch_currents(I, dz, alpha, Ra, Rb);
                      
            line(t/3600, Ia-Ib, 'LineStyle', linestyles{i}, 'LineWidth', 2, ...
                             'Color', 'k', ...
                             'DisplayName', sprintf('$R_a$ = %g m$\\Omega$, $R_b$ = %g m$\\Omega$', Ra*1e3, Rb*1e3))


%           Use this to check the consistency between different equations
%             dI = dI_dynamics(t, I, Ra, Rb, Qa, Qb, alpha, dz0);
% 
%             line(t/3600, dI, 'LineStyle', linestyles{i}, 'LineWidth', 1, ...
%                              'Color', 'm', ...
%                              'DisplayName', sprintf('Analytic'))

        end

        ylabel('$\Delta I$ (A)', 'Interpreter', 'Latex')
        yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
        lh = legend('show'); set(lh, 'Location', 'NorthEast', 'Interpreter', 'Latex')
        box on;

        % Branch current dynamics
        nexttile(th, 2)
      
        ax2 = gca(); box on;

        for i = 1:numel(r_vec)
            
            Rb = r_vec(i)*Ra;

            [za, zb] = solve_z_dynamics(t, I, alpha, Ra, Rb, Qa, Qb, za0, zb0);
            dz = za - zb;
            [Ia, Ib] = solve_branch_currents(I, dz, alpha, Ra, Rb);

            line(t/3600, Ia, 'Color', 'r', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
            line(t/3600, Ib, 'Color', 'b', 'LineStyle', linestyles{i}, 'LineWidth', 2, 'Parent', ax2)
    
            Vta = U0 + alpha.*za - Ia.*Ra;
            
            % Reference Lines
            time_at_vmax = interp1(Vta, t/3600, 4.2);
            xline(time_at_vmax, 'Color', 'm', 'LineWidth', 1, 'LineStyle', linestyles{i}, ...
                'Parent', ax1, 'HandleVisibility', 'off')            
            xline(time_at_vmax, 'Color', 'm', 'LineWidth', 1, 'LineStyle', linestyles{i}, ...
                'Parent', ax2, 'HandleVisibility', 'off')

        end
        xlabel(ax2, 'Time (hrs)', 'Interpreter', 'Latex')
        ylabel(ax2, '$I$ (A)', 'Interpreter', 'Latex')
        legend(ax2, {'$I_a$', '$I_b$'}, 'Location', 'South')
        box on;

        text(time_at_vmax-0.5, 1.0, '$V_t = 4.2V$', 'Color', 'm', 'FontSize', 16, ...
              'Parent', ax1, 'HandleVisibility', 'off', 'Interpreter', 'Latex')

        saveas(fh, sprintf('figures/fig_I_dynamics_timeseries_%g.png', iqb))

    end

end