function fig_nonlin_map()
    % Build imbalance maps for the nonlinear case

    set_default_plot_settings()

    % Configure the simulation
    analysis_type       = 'nmc';
    xvar                = 'i';
    to_plot_detailed    = false;

    % Initial conditions
    za0 = 0.0;
    zb0 = 0.0;
    dz0 = za0 - zb0;

    switch analysis_type
        case 'nmc'
            ocv = load_ocv_fn('nmc');
            Vmax = 4.2;
        case 'lfp'
            ocv =  load_ocv_fn('lfp');
            Vmax = 3.6;
    end

    Ra   = 0.15 ;    % Ohms
    Qa   = 3 * 3600; % Amp-seconds

    % Simulation vectors
    if to_plot_detailed
        num_grid_points = 5;
    else
        num_grid_points = 31;
    end

    simulation_hours = 16;
    dt = 1;
    t = 0:dt:simulation_hours*3600; t = t';
    I = -1/1 * (Qa / 3600) * ones(size(t)); % Current
    I_cv = I/5;

    r_vec = linspace(1.5, 0.5, num_grid_points);
    q_vec = linspace(0.5, 1.5, num_grid_points);
    di_ss = zeros(numel(r_vec), numel(q_vec));
    dz_ss = zeros(numel(r_vec), numel(q_vec));
    di_max = zeros(numel(r_vec), numel(q_vec));
    dz_max = zeros(numel(r_vec), numel(q_vec));

    if to_plot_detailed
        fh = figure('Position', [500 100 900, 900]);
        th = tiledlayout(numel(r_vec), numel(q_vec), ...
            'Padding', 'none', 'TileSpacing', 'tight', ...
            'TileIndexing', 'rowmajor');
    
        tile_count = 0;
        axes = [];
    end

    % Loop over the parameter space
    for i = 1:numel(r_vec)
        for j = 1:numel(q_vec)

            tic
            
            fprintf('Running (r,q) = (%.2f, %.2f)...\n', r_vec(i), q_vec(j))
            
            Rb = Ra / r_vec(i);
            Qb = Qa / q_vec(j);

            p.Qa = Qa;
            p.Qb = Qb;
            p.Ra = Ra;
            p.Rb = Rb;

            res = run_discrete_time_simulation_complete(I, -I, ...
               I_cv, Qa, Qb, Ra, Rb, za0, zb0, ocv, ...
               3, Vmax);
    
            % Filter out unwanted states
            idx = res.state == 2;
            res.t(idx) = [];
            res.za(idx) = [];
            res.zb(idx) = [];
            res.Ia(idx) = [];
            res.Ib(idx) = [];
            res.Vt(idx) = [];

            assert(~isnan(res.Ia(end)))

            di_ss(i, j) = res.Ia(end) - res.Ib(end);
            dz_ss(i, j) = res.za(end) - res.zb(end);

            dI = abs(res.Ia - res.Ib);
            idx = find(dI == max(dI));
            dI_sign = sign(res.Ia(idx) - res.Ib(idx));

            dz = abs(res.za - res.zb);
            idx = find(dz == max(dz));
            dz_sign = sign(res.za(idx) - res.zb(idx));
              
            di_max(i, j) = dI_sign(1).*max(abs(res.Ia - res.Ib));
            dz_max(i, j) = dz_sign(1).*max(abs(res.za - res.zb));

            % Calculate imbalance bounds
            [condition, zbound_l2, zbound_linf, ibound_l2, ibound_linf] = ...
                solve_imbalance_bounds(ocv, p, t, I, dz0);

            % Check if condition is satisfied
            is_condition_satisfied = max(abs(I)) < condition;

            if to_plot_detailed

                tile_count = tile_count + 1;

                ax = nexttile(th, tile_count); box on; 

                if strcmpi(xvar, 'z')
                    xa = res.za; xb = res.zb;
                    legend_string = {'$z_1$', '$z_2$'};
                    label_string = '$z$';
                elseif strcmpi(xvar, 'i')
                    xa = res.Ia; xb = res.Ib;
                    legend_string = {'$I_1$', '$I_2$'};
                    label_string = '$I$(A)';
                end

                line(res.t/3600, xb, 'Color', 'b', 'LineWidth', 2, 'LineStyle', ':', 'DisplayName', '$I_1$')
                line(res.t/3600, xa, 'Color', 'r', 'LineWidth', 2, 'DisplayName', '$I_2$')

%                 line(res.t/3600, xa - xb, 'Color', 'r', 'LineWidth', 2, 'DisplayName', '$I_a$')
%                 yline(ibound_linf);
%                 yline(1, 'LineStyle', '-', 'Color', 'r')

                axes = [axes ; ax];
                
                if i == 1 && j == 1 ; legend(legend_string) ; end
%                     else; legend(sprintf('%g', is_condition_satisfied)); end
                if j == numel(q_vec) ; set(ax, 'YTickLabel', []) ; end
                if i == numel(r_vec) ; xlabel(sprintf('$Q_2/Q_1$=%.2g\n$t$ (hrs)', q_vec(j)), 'Interpreter', 'Latex')
                    else; set(ax, 'XTickLabel', []) ; end
                if j == 1; ylabel(sprintf('%s\n$R_2/R_1$=%.2g', label_string,r_vec(i)), 'Interpreter', 'Latex')
                    else; set(ax, 'YTickLabel', []) ; end
            end

            toc

        end

    end

    if to_plot_detailed
        linkaxes(axes, 'xy'); 
    end


    make_plot(q_vec, r_vec, dz_max, '$|\Delta z|_{\mathrm{max}}$')
    make_plot(q_vec, r_vec, di_max, '$|\Delta I|_{\mathrm{max}}$ (A)')
    
end

function make_plot(q_vec, r_vec, z_map, title_string)

    % Make the figure
    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    [C,h] = contour(q_vec, r_vec, z_map, 'LabelSpacing', 300, ...
                                    'ShowText', 'on', 'LineColor', 'k', ...
                                    'Parent', gca, ...
                                    'HandleVisibility', 'off');

    hold on;
    contour(q_vec, r_vec, z_map, [0 0], ...
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
    ylabel('$r$ = $R_2/R_1$', 'Interpreter', 'Latex')
    xlabel('$q$ = $Q_2/Q_1$', 'Interpreter', 'Latex')
    title(title_string, 'Interpreter', 'Latex')

%     saveas(fh, 'figures/fig_dz_ss_sensitivity_nonlin_lfp.png')
%     saveas(fh, 'figures/fig_dz_ss_sensitivity_nonlin_lfp.fig')

end


function add_marker_annotations()

    % Add marker annotations for specific points in the parameter space
    qq = [1, 0.8, 1, 0.7];
    rr = [1, 1.25, 1.5, 1.1];

    marker_vec = {'o', 's', '^', 'o'};

    for i = 4

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