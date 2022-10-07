function test_lifetime_sei_growth_models()
    % Demonstrate a few SEI growth models

    set_default_plot_settings_manuscript()
    % Smith2011

    k = 1e-3;
    t = linspace(0, 1, 1000);
    x = sqrt(2.*k.*t);

    figure();
    plot(t, x, 'DisplayName', 'Smith2011')
    xlabel('Time (years)')
    legend show
    ylabel('x')


    k1 = 1e-3;
    k2 = 2.*k1;
    t = linspace(0, 1, 1000);
    x1 = sqrt(2.*k1.*t);
    x2 = sqrt(2.*k2.*t);

    figure();
    plot(t, x1, t, x2)
    xlabel('Time (years)')
    legend('k1', 'k2')
    ylabel('x')

end