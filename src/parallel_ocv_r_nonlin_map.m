function parallel_ocv_r_nonlin_map()
    % Build steady-state imbalance maps for the nonlinear case
    %
    % This code will take some time to complete

    clear all; close all; clc

    set_default_plot_settings_manuscript()

    %% Delta I_ss plot

    % Initial conditions
    za0 = 0.85;
    zb0 = 0.80;
    ocv = load_ocv_fn('lfp');
       
    Vmax = 3.6;
    Qa   = 5 * 3600; % Amp-seconds
    Ra   = 0.05 ;    % Ohms

    % Simulation vectors
    simulation_hours = 3;
    dt = 1;
    t = 0:dt:simulation_hours*3600; t = t';
    I = -1/3 * (Qa / 3600) * ones(size(t)); % Current

    num_points = 20;
    r_vec = linspace(0.5, 1.5, num_points + 1);
    q_vec = linspace(0.5, 1.5, num_points);
    di_ss = zeros(numel(r_vec), numel(q_vec));
    dz_ss = zeros(numel(r_vec), numel(q_vec));

    % Loop over the parameter space
    for i = 1:numel(r_vec)
        for j = 1:numel(q_vec)
            
            fprintf('Running (r,q) = (%.2f, %.2f)...\n', r_vec(i), q_vec(j))
            
            Rb = Ra / r_vec(i);
            Qb = Qa / q_vec(j);

            res = run_discrete_time_simulation(t, I, Qa, Qb, Ra, Rb, ...
                                               za0, zb0, ocv, Vmax);

            assert(~isnan(res.Ia(end)))

            di_ss(i, j) = res.Ia(end) - res.Ib(end);
            dz_ss(i, j) = res.za(end) - res.zb(end);

        end
    end
   
    % Make the figure
    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    [C,h] = contour(q_vec, r_vec, di_ss, 'LabelSpacing', 300, ...
                                    'ShowText', 'on', 'LineColor', 'k', ...
                                    'Parent', gca, ...
                                    'HandleVisibility', 'off');

    hold on;
    contour(q_vec, r_vec, di_ss, [0 0], ...
            'LineWidth', 5, 'Color', 'k', ...
            'LineStyle', '--', 'HandleVisibility', 'off')

    clabel(C, h, 'FontSize', 24, ...
                 'Color', 'k', ...
                 'FontName', 'Times New Roman', ...
                 'Interpreter', 'Latex')

    set(findobj(gca, 'Type', 'patch', 'UserData', 2), ...
        'EdgeColor', [0 0 0])

    box on;

    % Add marker annotations for specific points in the parameter space
    add_marker_annotations()


    lh = legend('show'); set(lh, 'Interpreter', 'Latex', ...
                                 'Location', 'SouthWest')

    yline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xlim([min(r_vec), max(r_vec)])
    ylim([min(q_vec), max(q_vec)])
    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\Delta I_{f}$ (A)', 'Interpreter', 'Latex')

    saveas(fh, 'figures/fig_di_ss_sensitivity_nonlin_lfp.fig')
    saveas(fh, 'figures/fig_di_ss_sensitivity_nonlin_lfp.png')

    % Make the figure
    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    [C,h] = contour(q_vec, r_vec, dz_ss, ...
                                    'LabelSpacing', 300, ...
                                    'ShowText', 'on', 'LineColor', 'k', ...
                                    'Parent', gca, ...
                                    'HandleVisibility', 'off');

    hold on;
    contour(q_vec, r_vec, dz_ss, [0 0], 'LineWidth', 5, 'Color', 'k', ...
                'LineStyle', '--', 'HandleVisibility', 'off')

    clabel(C, h, 'FontSize', 24, ...
                 'Color', 'k', ...
                 'FontName', 'Times New Roman', ...
                 'Interpreter', 'Latex')

    set(findobj(gca, 'Type', 'patch', 'UserData', 2), ...
        'EdgeColor', [0 0 0])

    box on;

    add_marker_annotations()


    lh = legend('show'); set(lh, 'Interpreter', 'Latex', ...
                                 'Location', 'SouthWest')

    yline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xline(1, 'LineStyle', ':', 'HandleVisibility', 'off')
    xlim([min(r_vec), max(r_vec)])
    ylim([min(q_vec), max(q_vec)])
    ylabel('$r$ = $R_a/R_b$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_a/Q_b$', 'Interpreter', 'Latex')
    title('$\Delta z_{f}$', 'Interpreter', 'Latex')

    saveas(fh, 'figures/fig_dz_ss_sensitivity_nonlin_lfp.png')
    saveas(fh, 'figures/fig_dz_ss_sensitivity_nonlin_lfp.fig')


    keyboard

    
end

function add_marker_annotations()

    % Add marker annotations for specific points in the parameter space
    qq = [1, 0.8, 1];
    rr = [1, 1.25, 1.5];

    marker_vec = {'o', 's', '^'};

    for i = 2

        q = qq(i);
        r = rr(i);

        line(q, r, ...
        'Marker', marker_vec{i},  ...
        'Color', [0 0.5 0], ...
        'MarkerFaceColor', [0 0.5 0], ...
        'LineStyle', 'none', ...
        'MarkerSize', 12, ...
        'DisplayName', sprintf('$(q,r)=(%g,%g)$', q, r));

    end

end