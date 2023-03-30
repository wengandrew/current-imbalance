function fig_ocv_functions()
    % Build reference curves highlighting the OCV functions
    set_default_plot_settings()

       
    nmc = load_ocv_fn('nmc');
    lfp = load_ocv_fn('lfp');
    
    nmc_affine = load_ocv_fn('nmc-affine');
    lfp_affine = load_ocv_fn('lfp-affine');

    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    xx = linspace(0, 1, 1000);

    line(xx, nmc(xx), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'b', 'DisplayName', 'NMC/Gr, nonlinear')
    line(xx, lfp(xx), 'LineWidth', 2, 'LineStyle', '-', 'Color', 'r', 'DisplayName', 'LFP/Gr, nonlinear')
    line(xx, nmc_affine(xx), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b', 'DisplayName', 'NMC/Gr, affine')
    line(xx, lfp_affine(xx), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r', 'DisplayName', 'LFP/Gr, affine')
    xlabel('$z$', 'Interpreter', 'Latex')
    ylabel('$U(z)$', 'Interpreter', 'Latex')
    box on
    lh = legend('show'); set(lh, 'Location', 'NorthWest')

end