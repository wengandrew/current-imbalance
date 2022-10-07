function parallel_ocv_r_nonlin_comparison()
    % Workspace for running discrete-time simulations of the parallel OCV-R
    % model. 
    %
    % This work serves to verify that the analytic equations for current
    % and SOC imbalance are in fact correct. The explicit numerical
    % simulation will also help explore what happens if we have non-linear
    % OCV functions or non-linear R functions.

    U0 = 3.0;
    Vmax = 4.2;
    alpha = 1.2;

    % Declare OCV functions
    ocv_lin = @(z) U0 + alpha * z;
    ocv_lfp = load_ocv_fn('lfp');
    ocv_nmc = load_ocv_fn('nmc');
    
    ocv_nmc_piecewise = @(z) ocv_lin_piecewise(z);

    za0 = 0.85;
    zb0 = 0.80;

    Qa = 5 * 3600; % Amp-seconds
    Ra = 0.05 ; % Ohms

    q = 1;
    r = 1.5;

    Rb = Ra/r;
    Qb = Qa/q;

    simulation_hours = 3;

    dt = 1;
    torig = 0:dt:simulation_hours*3600; torig = torig';
    I = -1/3 * (Qa / 3600) * ones(size(torig)); % Current

    nmcpw = run_discrete_time_simulation(torig, I, Qa, Qb, Ra, Rb, ...
        za0, zb0, ocv_nmc_piecewise, Vmax);
    
    % Analytical results from the paper to verify against
    [zaa, zbb] = solve_z_dynamics(torig, I, alpha, Ra, Rb, Qa, Qb, za0, zb0);
    dzz = zaa - zbb;
    [Iaa, Ibb] = solve_branch_currents(I, dzz, alpha, Ra, Rb);
    dII = Iaa - Ibb;
    Vtt = U0 + alpha*zaa - Iaa * Ra;

    % Truncate the data
    dzz(nmcpw.idx) = [];
    dII(nmcpw.idx) = [];
    zaa(nmcpw.idx) = [];
    zbb(nmcpw.idx) = [];
    Iaa(nmcpw.idx) = [];
    Ibb(nmcpw.idx) = [];
    Vtt(nmcpw.idx) = [];

    % Visualize
    set_default_plot_settings_manuscript()

    fh = figure('Position', [500 100 800*1.5 500*1.5]);
    th = tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'tight'); 
   
     
    ax_dz = nexttile(th, 1); box on;
    ax_dz = gca(); % set(ax_dz,'XTickLabel',[])
    ylabel(ax_dz, '$\Delta z$', 'Interpreter', 'Latex')
    xlabel(ax_dz, 'Time (hrs)', 'Interpreter', 'Latex')

    ax_z = nexttile(th, 4); box on;
    ax_z = gca();
    ylabel(ax_z, '$z$', 'Interpreter', 'Latex')    
    xlabel(ax_z, 'Time (hrs)')

    ax_di = nexttile(th, 2); box on;
    ylabel(ax_di, '$\Delta I$ (A)', 'Interpreter', 'Latex')
    xlabel(ax_di, 'Time (hrs)')

    ax_i = nexttile(th, 5); box on;
    ylabel(ax_i, '$I$ (A)', 'Interpreter', 'Latex')
    xlabel(ax_i, 'Time (hrs)')

    ax_vt = nexttile(th, 6); box on;
    ylabel(ax_vt, '$V_t$ (V)', 'Interpreter', 'Latex')
    xlabel(ax_vt, 'Time (hrs)')

    ax_voc = nexttile(th, 3); box on;
    ylabel(ax_voc, '$V_{OC} (V)$', 'Interpreter', 'Latex')
    xlabel(ax_voc, '$z$', 'Interpreter', 'Latex')

    % Plot the OCV curves
    zz = linspace(0, 1, 1000);
    line(zz, ocv_lfp(zz), 'Color', 'b', 'Parent', ax_voc)
    line(zz, ocv_lin(zz), 'Color', 'k', 'LineStyle', '-', ...
                'DisplayName', 'Analytic', 'Parent', ax_voc)
    line(zz, ocv_nmc(zz), 'Color', 'r', ...
                'DisplayName', 'NMC', 'Parent', ax_voc)
    line(zz, ocv_nmc_piecewise(zz), 'Color', 'r', ...
                'DisplayName', 'NMCLinPiecewise', 'LineStyle', '--', 'Parent', ax_voc)
    legend(ax_voc, 'LFP', 'NMC', 'Linear')
    legend(ax_voc, 'show')

    line(nmcpw.t./3600, dzz, 'Color', 'k', ...
        'DisplayName', 'Analytic', 'LineStyle', '-', 'Parent', ax_dz)

    line(nmcpw.t./3600, nmcpw.za, 'Color', 'r', 'LineStyle', '--', 'Parent', ax_z)
    line(nmcpw.t./3600, nmcpw.zb, 'Color', 'r', 'LineStyle', '--', 'Parent', ax_z)
    line(nmcpw.t./3600, zaa, 'Color', 'k', 'LineStyle', '-', 'Parent', ax_z)
    line(nmcpw.t./3600, zbb, 'Color', 'k', 'LineStyle', '-', 'Parent', ax_z)

    line(nmcpw.t./3600, nmcpw.Ia - nmcpw.Ib, 'Color', 'r', 'LineStyle', '--', 'Parent', ax_di)
    line(nmcpw.t./3600, dII, 'Color', 'k', 'LineStyle', '-', 'Parent', ax_di)

    line(nmcpw.t./3600, nmcpw.Ia, 'Color', 'r', 'LineStyle', '--', 'Parent', ax_i)
    line(nmcpw.t./3600, nmcpw.Ib, 'Color', 'r', 'LineStyle', '--', 'Parent', ax_i)
    line(nmcpw.t./3600, Iaa, 'Color', 'k', 'LineStyle', '-', 'Parent', ax_i)
    line(nmcpw.t./3600, Ibb, 'Color', 'k', 'LineStyle', '-', 'Parent', ax_i)

    line(nmcpw.t./3600, nmcpw.Vt, 'Color', 'r', 'LineStyle', '--', 'Parent', ax_vt)
    line(nmcpw.t./3600, Vtt, 'Color', 'k', 'LineStyle', '-', 'Parent', ax_vt)


    % Do the LFP scenario
    Vmax = 3.6;
    lfp = run_discrete_time_simulation(torig, I, Qa, Qb, Ra, Rb, ...
        za0, zb0, ocv_lfp, Vmax);
    line(lfp.t./3600, lfp.za - lfp.zb, 'Color', 'b', 'DisplayName', 'LFP', 'Parent', ax_dz)
    line(lfp.t./3600, lfp.za, 'Color', 'b', 'Parent', ax_z)
    line(lfp.t./3600, lfp.zb, 'Color', 'b', 'Parent', ax_z)
    line(lfp.t./3600, lfp.Ia - lfp.Ib, 'Color', 'b', 'Parent', ax_di)
    line(lfp.t./3600, lfp.Ia, 'Color', 'b', 'Parent', ax_i)
    line(lfp.t./3600, lfp.Ib, 'Color', 'b', 'Parent', ax_i)
    line(lfp.t./3600, lfp.Vt, 'Color', 'b', 'Parent', ax_vt)

    % Do the NMC scenario
    Vmax = 4.2;
    nmc = run_discrete_time_simulation(torig, I, Qa, Qb, Ra, Rb, ...
        za0, zb0, ocv_nmc, Vmax);

    line(nmc.t./3600, nmc.za - nmc.zb, 'Color', 'r', 'DisplayName', 'NMC', 'Parent', ax_dz)
    line(nmc.t./3600, nmc.za, 'Color', 'r', 'Parent', ax_z)
    line(nmc.t./3600, nmc.zb, 'Color', 'r', 'Parent', ax_z)
    line(nmc.t./3600, nmc.Ia - nmc.Ib, 'Color', 'r', 'Parent', ax_di)
    line(nmc.t./3600, nmc.Ia, 'Color', 'r', 'Parent', ax_i)
    line(nmc.t./3600, nmc.Ib, 'Color', 'r', 'Parent', ax_i)
    line(nmc.t./3600, nmc.Vt, 'Color', 'r', 'Parent', ax_vt)


    line(nmcpw.t./3600, nmcpw.za - nmcpw.zb, 'Color', 'r', ...
        'DisplayName', 'NMCLinPiecewise', 'LineStyle', '--', 'Parent', ax_dz)

    lh = legend(ax_dz, 'show'); set(lh, 'Location', 'NorthEast')
    legend(ax_dz, 'LinPiecewise', 'Analytic', 'NMC')

    saveas(fh, 'figures/fig_nonlinear_study.png')

    keyboard

end

function ocv = ocv_lin_piecewise(z_vec)
    % Piece-wise affine OCV function for an NMC/Graphite system.
    % 
    % The equations are assembled manually by picking 4 points 
    % from an NMC/grahpite curve and defining affine functions
    % connecting them together.

    [coeffs1, ~] = polyfit([0, 0.14014], [2.75973 3.43615], 1);
    [coeffs2, ~] = polyfit([0.140141, 0.834835], [3.43615, 4.07314], 1);
    [coeffs3, ~] = polyfit([0.834835, 0.935936], [4.07314, 4.11217], 1);
    [coeffs4, ~] = polyfit([0.935936, 1], [4.11217, 4.19999], 1);

    for i = 1:numel(z_vec)

        z = z_vec(i);

        if (z >= 0) && (z < 0.14014)
            alpha = coeffs1(1);
            U0 = coeffs1(2);
        elseif (z >= 0.14014) && (z < 0.834835)
            alpha = coeffs2(1);
            U0 = coeffs2(2);
        elseif (z >= 0.834835) && (z < 0.935936)
            alpha = coeffs3(1);
            U0 = coeffs3(2);
        else
            alpha = coeffs4(1);
            U0 = coeffs4(2);
        end

        ocv(i) = U0 + alpha * z;

    end

end