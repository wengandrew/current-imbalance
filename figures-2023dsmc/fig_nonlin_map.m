function fig_nonlin_map()
    % Build imbalance maps for the nonlinear case

    set_default_plot_settings()

    % Configure the simulation
    analysis_type       = 'lfp';
    xvar                = 'z'; % z (soc) or i (current)
    to_plot_detailed    = false;
    to_plot_imbalance   = true;

    % Initial conditions
    za0 = 0.0;
    zb0 = 0.0;
    dz0 = za0 - zb0;

    switch analysis_type
        case 'nmc'
            ocv = load_ocv_fn('nmc');
            Vmin = 3.0;
            Vmax = 4.2;
            affine_name = 'nmc-affine';
        case 'lfp'
            ocv =  load_ocv_fn('lfp');
            Vmin = 3.0;
            Vmax = 3.6;
            affine_name = 'lfp-affine';
    end

    Ra   = 0.05 ;    % Ohms
    Qa   = 3 * 3600; % Amp-seconds

    % Simulation vectors
    if to_plot_detailed
        num_grid_points = 4;
    else
        num_grid_points = 60;
    end

    ocv_lin = load_ocv_fn(affine_name);
    Vmax_affine = ocv_lin(1);
    Vmin_affine = ocv_lin(0);
    alpha = Vmax_affine - Vmin_affine;

    simulation_hours = 160;
    dt = 1;
    t = 0:dt:simulation_hours*3600; t = t';
    I = -1/1 * (Qa / 3600) * ones(size(t)); % Current
    I_cv = I(1)/5;

    r_vec = linspace(1.5, 0.5, num_grid_points);
    q_vec = linspace(0.5, 1.5, num_grid_points);
    di_ss = zeros(numel(r_vec), numel(q_vec));
    dz_ss = zeros(numel(r_vec), numel(q_vec));
    di_max = zeros(numel(r_vec), numel(q_vec));
    dz_max = zeros(numel(r_vec), numel(q_vec));

    di_max_sign = zeros(numel(r_vec), numel(q_vec));
    dz_max_sign = zeros(numel(r_vec), numel(q_vec));

    if to_plot_detailed
        fh = figure('Position', [500 100 1200, 700]);
        th = tiledlayout(numel(r_vec), numel(q_vec), ...
            'Padding', 'none', 'TileSpacing', 'none', ...
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

            % Nonlinear solution
            res = run_discrete_time_simulation_multicycle(I, -I, ...
               I_cv, Qa, Qb, Ra, Rb, za0, zb0, ocv, ...
               3, Vmax);

            kappa = 1 / alpha * (Ra*Qa - Rb*Qb) / (Qa + Qb);


            % Affine solution
            res_aff = solve_z_dynamics_cccv_complete(t, I, -I, ...
            I_cv, Ra, Rb, Qa, Qb, za0, zb0, ocv_lin, Vmin, Vmax);

%             % Filter out unwanted states
%             idx = find(res_aff.t >= res_aff.t_chg_cv);
%             res_aff.t(idx) = [];
%             res_aff.za(idx) = [];
%             res_aff.zb(idx) = [];
%             res_aff.Ia(idx) = [];
%             res_aff.Ib(idx) = [];
%             res_aff.Vt(idx) = [];
%     
%             % Filter out unwanted states
%             idx = res.state == 2;
%             res.t(idx) = [];
%             res.za(idx) = [];
%             res.zb(idx) = [];
%             res.Ia(idx) = [];
%             res.Ib(idx) = [];
%             res.Vt(idx) = [];

            assert(~isnan(res.Ia(end)))

            di_ss(i, j) = res.Ia(end) - res.Ib(end);
            dz_ss(i, j) = res.za(end) - res.zb(end);

            dI = abs(res.Ia - res.Ib);
            idx = find(dI == max(dI));
            dI_sign = sign(res.Ia(idx) - res.Ib(idx));

            dz = abs(res.za - res.zb);
            idx = find(dz == max(dz));
            dz_sign = sign(res.za(idx) - res.zb(idx));

            di_max(i, j) = max(abs(res.Ia - res.Ib));
            dz_max(i, j) = max(abs(res.za - res.zb));

            di_max_sign(i, j) = dI_sign(end);
            dz_max_sign(i, j) = dz_sign(end);

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
                    xaa = res_aff.za; xbb = res_aff.zb;
                    labels = {'$z_1$', '$z_2$'};
                    label_string = '$z$';
                    imbalance_string = '$\Delta z$';
                    boundvar = zbound_linf;
                    if boundvar > 1; boundvar = 1; end
                    boundvar_aff = kappa*I(1);
                elseif strcmpi(xvar, 'i')
                    xa = res.Ia; xb = res.Ib;
                    xaa = res_aff.Ia; xbb = res_aff.Ib;
                    labels = {'$I_1$', '$I_2$'};
                    label_string = '$I$(A)';
                    imbalance_string = '$\Delta I$';
                    boundvar = ibound_linf;
                    boundvar_aff = (Qa - Qb) / (Qa + Qb) * I(1);
                end


                % Define time offsets so that the affine and the non-linear simulations
                % will share common points at the start and end of the CV hold phase.
                t0 = res_aff.t_chg_cc;
                t1 = res_aff.t_chg_cv;
                t3 = res.t(find(res.state==1, 1, 'first'));
                t2 = res.t(end);
                idxa = res.t < t2;
            
                xline(res_aff.t_chg_cc./3600 - (t0-t3)/3600, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')
    
                if to_plot_imbalance
                    labels = {imbalance_string};
                    label_string = imbalance_string;
                    line(res_aff.t/3600 - (t0-t3)/3600, xaa - xbb, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', [imbalance_string ' (Affine)'])
                    yline(boundvar_aff, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Bound, Affine ($\kappa$I)')
                    line(res.t/3600, xa - xb, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', [imbalance_string ' (NMC/Gr)'])
                    yline(boundvar, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'Bound, Nonlinear')

                else
                    line(res_aff.t/3600 - (t0-t3)/3600, xbb, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', [labels{1} ' (affine)'])
                    line(res_aff.t/3600 - (t0-t3)/3600, xaa, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', [labels{2} ' (affine)'])
                    line(res.t/3600, xb, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', [labels{1} ' (NMC/Gr)'])
                    line(res.t/3600, xa, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', [labels{2} ' (NMC/Gr)'])
                    yline(0, 'LineStyle', ':', 'Color', 'k', 'LineWidth', 1, 'HandleVisibility', 'off')

                end

%                 ylim([-1 1])
                axes = [axes ; ax];
                
                if i == 1 && j == 1 ; lh = legend('show') ; end
%                     else; legend(sprintf('%g', is_condition_satisfied)); end
                if j == numel(q_vec) 
                    set(ax, 'YTickLabel', []) ; 
                end
                if i == numel(r_vec) ; xlabel(sprintf('$t$ (hrs) \n $Q_2/Q_1$=%.2g', q_vec(j)), 'Interpreter', 'Latex')
                    else; set(ax, 'XTickLabel', []) ; end
                if j == 1; ylabel(sprintf('$R_2/R_1$=%.2g\n %s', r_vec(i),label_string), 'Interpreter', 'Latex')
                    else
                        set(ax, 'YTickLabel', []) ; 
                end
            end

            toc

        end

    end

    if to_plot_detailed
        linkaxes(axes, 'xy'); 
    end


    make_plot(q_vec, r_vec, dz_max, '$|\Delta z|_{\mathrm{max}}$')
    make_plot(q_vec, r_vec, di_max, '$|\Delta I|_{\mathrm{max}}$ (A)')
keyboard
    
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
    qq = [1, 0.8, 1, 0.7, 0.7];
    rr = [1, 1.25, 1.5, 1.1, 1/0.7];

    marker_vec = {'o', 's', '^', 'o', 's'};

    for i = 4:5

        q = qq(i);
        r = rr(i);

        line(q, r, ...
        'Marker', marker_vec{i},  ...
        'Color', [0 0.5 0], ...
        'MarkerFaceColor', [0 0.5 0], ...
        'LineStyle', 'none', ...
        'MarkerSize', 12, ...
        'DisplayName', sprintf('$(q,r)=(%.2f,%.2f)$', q, r));

    end

end