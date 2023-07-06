function aw20230303_model_vs_measured(plottype)
    % Run a full charge CCCV and discharge CC simulation
    %
    % Parameters
    % ----------
    % plottype: 'c10' or 'c4' only

    set_default_plot_settings()

    tbla = readtable('data/params_cell151805.csv');
    tblb = readtable('data/params_cell152098.csv');

    ocv_a = griddedInterpolant(tbla.soc/100, tbla.ocv, 'linear', 'linear');
    ocv_b = griddedInterpolant(tblb.soc/100, tblb.ocv, 'linear', 'linear');

    % Initialize model parameters
    
    % Add general support for SOC-dependent resistances.
    % But to make the results consistent with the proposed OCV-R model, 
    % we will actually make these look-up functions return static values
    % for DCR which equal the SOC-averaged DCR value for each cell.
    Ra = griddedInterpolant(tbla.soc/100, ones(size(tbla.dcr))*mean(tbla.dcr), 'linear', 'linear');
    Rb = griddedInterpolant(tblb.soc/100, ones(size(tblb.dcr))*mean(tblb.dcr), 'linear', 'linear');
    
    % Uncomment out the lines below to get the non-linear resistance table
    % lookups.
    %     Ra = griddedInterpolant(tbla.soc/100, tbla.dcr, 'linear', 'linear');
    %     Rb = griddedInterpolant(tblb.soc/100, tblb.dcr, 'linear', 'linear');

    Qa = 1.83 * 3600; % As - measured from C/20 discharge capacity
    Qb = 1.93 * 3600; % As - measured from C/20 discharge capacity
    
   
    za0 = 0.00325; % Tune to start at around 3.035V, to match the data
    zb0 = 0.00325; % Tune to start at around 3.035V, to match the data

    Vmax = 4.2;
    Vmin = 3.0;
    U0 = 3.0;
    
    alpha = Vmax - Vmin;

    % Get the raw data in there
    if strcmpi(plottype, 'c10')
        data = readtable('data/expt_c_over_10.csv');
        data.test_time_s = data.test_time_s - 0.98294*3600 + 0.146*3600;
        hours = 10;
    elseif strcmpi(plottype, 'c4')
        data = readtable('data/expt_c_over_4.csv');
        data.test_time_s = data.test_time_s - 0.98294*3600;
        hours = 4;
    end

    current_target = - 2 * 2.5 / hours; % Total current / Amperes

    %% Test the analytic solution
    ocv_lin = @(z) U0 + alpha * z;

    % Initialize simulation parameters
    t = linspace(0, (hours+5)*3600, 1.0e6)';
    I_chg = +current_target*ones(size(t)); % applied current (A) 
    I_dch = -current_target*ones(size(t));
    I_cutoff = 2.5/30;
    
    res_lsim = solve_z_dynamics_cccv_complete(t, I_chg, I_dch, ...
        I_cutoff, mean(Ra(0:0.01:1)), mean(Rb(0:0.01:1)), Qa, Qb, za0, zb0, ocv_lin);

    res_disc = run_discrete_time_simulation_complete_nonlinr(I_chg, I_dch, ...
        I_cutoff, Qa, Qb, Ra, Rb, za0, zb0, ocv_a, ocv_b, Vmin, Vmax);


    % Current imbalance
    fh = figure('Position', [500 100 650 600]);
    th = tiledlayout(2, 1, 'Padding', 'none', 'TileSpacing', 'compact'); 
    
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);

    tt = data.test_time_s;
    i1 = sgolayfilt(-data.current_1, 3, 21);
    i2 = sgolayfilt(-data.current_2, 3, 21);

    line(tt./3600, i2, 'Color', 'b', 'Parent', ax1, ...
        'LineWidth', 2, 'DisplayName', 'Experiment, Cell 1')    
    line(tt./3600, i1, 'Color', 'r', 'Parent', ax1, ...
        'LineWidth', 2, 'DisplayName', 'Experiment, Cell 2')
    line(res_disc.t./3600, res_disc.Ib, 'Color', 'b', ...
        'LineWidth', 2, 'LineStyle', '--', 'Parent', ax1, ...
        'DisplayName', 'Model, Cell 1')
    line(res_disc.t./3600, res_disc.Ia, 'Color', 'r', ...
        'LineWidth', 2, 'LineStyle', '--', 'Parent', ax1, ...
        'DisplayName', 'Model, Cell 2')

    xlabel(ax1, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax1, '$I$ (A)', 'Interpreter', 'Latex')
    lh = legend(ax1, 'show'); set(lh, 'Location', 'northeast')
    
    % Voltage
    ax2 = nexttile(th, 2); box on; 

    line(data.test_time_s./3600, data.voltage_v, 'LineWidth', 2, ...
        'Color', 'k', 'Parent', ax2, 'DisplayName', 'Experiment, $V_t$')
    line(res_disc.t./3600, res_disc.Vt, 'LineWidth', 2, 'Color', 'k', ...
        'LineStyle', '--', 'Parent', ax2, 'DisplayName', 'Model, $V_t$')
    line(res_disc.t./3600, ocv_b(res_disc.zb), 'LineWidth', 2, ...
        'Color', 'b', 'LineStyle', '--', 'Parent', ax2, 'DisplayName', 'Model, $U_1$')
    line(res_disc.t./3600, ocv_a(res_disc.za), 'LineWidth', 2, ...
        'Color', 'r', 'LineStyle', '--', 'Parent', ax2, 'DisplayName', 'Model, $U_2$')
    ylim([3 4.25])


    xlabel(ax2, '$t$ (hrs)', 'Interpreter', 'Latex')
    ylabel(ax2, '$V$ (V)', 'Interpreter', 'Latex')
    lh = legend(ax2, 'show'); set(lh, 'Location', 'northeast')

    linkaxes([ax1, ax2], 'x')

    if strcmpi(plottype, 'c10')
        xlim(ax2, [0, 15.47])
        line([0 15.47], [0 0], 'Color', [0.5, 0.5, 0.5], ...
            'LineStyle', ':', 'HandleVisibility', 'off', 'Parent', ax1)

    elseif strcmpi(plottype, 'c4')
        xlim(ax2, [0, 7.39])
        line([0 7.39], [0 0], 'Color', [0.5, 0.5, 0.5], ...
            'LineStyle', ':', 'HandleVisibility', 'off', 'Parent', ax1)
    end

end
