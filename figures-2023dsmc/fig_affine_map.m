function fig_affine_map()
    %% Andrew Weng
    
    set_default_plot_settings()

    % Initialize model parameters
    Ra = 0.150; % ohms
    Qa = 3 * 3600; % As
    alpha = 1.2;
    za0 = 0.0;
    zb0 = 0.0;
    U0  = 3.0;
    current_target = -Qa / (1 * 3600);

    % Initialize simulation vectors
    t = linspace(0, 2.5*3600, 1e4)';
    I = current_target*ones(size(t)); % applied current (A) 

    %% Make the plots
    plot_dz_ss_sensitivity(Qa, Ra, current_target, alpha)
    plot_di_ss_sensitivity(Qa, Ra, current_target)
    
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

    ylabel('$r$ = $R_2/R_1$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_2/Q_1$', 'Interpreter', 'Latex')
    title('$\tau$ (hrs)', 'Interpreter', 'Latex')
%     saveas(fh, 'figures/fig_time_constant.png')

end

function plot_dz_ss_sensitivity(Qa, Ra, current_target, alpha)

    set_default_plot_settings()

    rr = linspace(0.5, 1.5, 101);
    qq = linspace(0.5, 1.5, 100);

    % Note: this can be vectorized
    for i = 1:numel(rr)
        for j = 1:numel(qq)
            Rb = Ra / rr(i);
            dz_ss(i, j) = Rb/alpha * (rr(i)*qq(j) - 1) / (qq(j) + 1) * current_target;
        end
    end

    c1 = [0.7, 0.7, 0.7];
    c2 = [1, 1, 1];
    c3 = [0.4, 0.4, 0.4];

    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    grid_points = [-0.3:0.05:0.3];

    [C,h] = contourf(qq, rr, dz_ss, grid_points, ...
                                 'ShowText', 'on', 'LineColor', c3, ...
                                 'Parent', gca, 'LabelSpacing', 270, ...
                                 'HandleVisibility', 'off');
    
    clabel(C, h, grid_points, 'FontSize', 15, ...
                                      'Color', c3, 'Interpreter', 'Latex')

    colormap([c2; c2; c2; c2; c1; c1; c1; c1; c1; c1; c1; c1; c1; c1])
    
    line([0.5 1.5], [1 1], 'Color', [0.3 0.3 0.3], ...
                           'LineStyle', ':', ...
                           'LineWidth', 1.0, 'HandleVisibility', 'off')

    line([1 1], [0.5 1.5], 'Color', [0.3 0.3 0.3], ...
                           'LineStyle', ':', ...
                           'LineWidth', 1.0, 'HandleVisibility', 'off')

    rectangle('Position', [0.5 1 0.5 0.5], 'EdgeColor',' r', ...
              'LineWidth', 3, 'LineStyle', '-', 'FaceColor', [1 0 0 0.07])

    line(rr, 1./rr, ...
        'LineWidth', 3, ...
        'Color',' k', ...
        'LineStyle', '--', ...
        'Parent', gca, ...
        'DisplayName', '$Q_2R_2 = Q_1R_1$'); hold all;

    text(1.25, 0.86, '$Q_2R_2 = Q_1R_1$', ...
                     'Interpreter', 'Latex', 'FontSize', 22', 'Rotation', -33)
    text(0.60, 0.90, '$z_{2} > z_{1}$', ...
                     'Interpreter', 'Latex', 'FontSize', 22, 'BackgroundColor', c1)
    text(1.10, 1.10, '$z_{2} < z_{1}$', ...
                     'Interpreter', 'Latex', 'FontSize', 22, 'BackgroundColor', c2)

    q_vec = [0.7, 0.8, 1];
    r_vec = [1.1, 1.25, 1.25];
    marker_vec = {'o', 's', '^'};

    for i = 1

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

    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])

    xlim([min(rr), max(rr)])
    ylim([min(qq), max(qq)])

    ylabel('$R_2/R_1$', 'Interpreter', 'Latex')
    xlabel('$Q_2/Q_1$', 'Interpreter', 'Latex')
    title('$\Delta z_{\mathrm{ss}}$', 'Interpreter', 'Latex')
    box on;

%     saveas(fh, 'figures/fig_dz_ss_sensitivity.png')

end

function plot_di_ss_sensitivity(Qa, Ra, current_target)

    set_default_plot_settings()

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

    c1 = [0.7, 0.7, 0.7];
    c2 = [1, 1, 1];
    c3 = [0.4, 0.4, 0.4];

    [C,h] = contourf(qq, rr, di_ss, [-1:0.2:1], ...
                                 'ShowText', 'off', 'LineColor', c3, ...
                                 'Parent', gca, 'LabelSpacing', 250, ...
                                 'HandleVisibility', 'off');

    
    clabel(C, h, [-1:0.2:1], 'FontSize', 15, 'Color', c3, 'Interpreter', 'Latex')

    colormap([c1; c1; c2; c2; c2; c2; c2])

    line([0.5 1.5], [1 1], 'Color', [0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth', 1.0, 'HandleVisibility', 'off')
    line([1 1], [0.5 1.5], 'Color', [0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth', 1.0, 'HandleVisibility', 'off')

    rectangle('Position', [0.5 1 0.5 0.5], 'EdgeColor',' r', ...
              'LineWidth', 3, 'LineStyle', '-', 'FaceColor', [1 0 0 0.07])


    line(ones(size(qq)), qq, ...
    'LineWidth', 3, ...
    'Color',' k', ...
    'LineStyle', '--', ...
    'Parent', gca, ...
    'DisplayName', '$q=1$'); hold all;

    q_vec = [0.7, 0.8, 1];
    r_vec = [1.1, 1.25, 1.25];
    marker_vec = {'o', 's', '^'};

    text(0.57, 0.90, '$|I_{ss,2}| < |I_{ss,1}|$', 'Interpreter', 'Latex', 'FontSize', 22, 'BackgroundColor', c2)
    text(1.07, 1.10, '$|I_{ss,2}| > |I_{ss,1}|$', 'Interpreter', 'Latex', 'FontSize', 22, 'BackgroundColor', c1)

    for i = 1

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

    set(findobj(gca,'Type','patch','UserData',2),'EdgeColor',[0 0 0])

    xlim([min(rr), max(rr)])
    ylim([min(qq), max(qq)])

    ylabel('$R_2/R_1$', 'Interpreter', 'Latex')
    xlabel('$Q_2/Q_1$', 'Interpreter', 'Latex')
    title('$\Delta I_{\mathrm{ss}}$ (A)', 'Interpreter', 'Latex')
    box on;

%     saveas(fh, 'figures/fig_di_ss_sensitivity.png')

end