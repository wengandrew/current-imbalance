function fig_sei_derivation()
    % A figure to explain how the SEI incremental loss update equation is
    % derived

    set_default_plot_settings();

    xx = [0, 1, 2, 3];

    k = 0.05;
    p = 0.5;

    fh = figure('Position', [500 100 400*1.5 400*1.5]);

    L = @(x) k .* xx .^ p;
    
    elem(xx, L(xx), 'lw', 2, 'c', 'k', 'm', 'o', 'ms', 10)
    xlabel('$k$', 'Interpreter', 'LaTeX')
    ylabel('$L_k$', 'Interpreter', 'LaTeX')

end