function parallel_ocv_r_galvanostatic()
    % Simple test for galvanostatic mode (CV) simulations with the parallel OCV-R model

    Uf = 4.2;
    U0 = 3.0;
    alpha = 1.2;

    Qa = 5.00 * 3600; % A-s
    Ra = 0.05;
    Rb = 0.07;

    r = 1.2;
    q = 0.8;

    Rb = Ra/r;
    Qb = Qa/q;

    za0 = 0.9;
    zb0 = 0.94;

    A = [-alpha/(Qa*Ra) 0 ; ...
         0 -alpha/(Qb*Rb)];

    B = [1/(Qa*Ra) ; 1/(Qb*Rb)];

    C = eye(2);

    D = 0;

    t = 0:0.1:5*3600;
    u = ones(size(t)) * (Uf - U0);

    sys = ss(A, B, C, D);
    out = lsim(sys, u, t, [za0;zb0]);

    za = out(:, 1);
    zb = out(:, 2);

    Ia = ( alpha*za - (Uf - U0) ) / Ra;
    Ib = ( alpha*zb - (Uf - U0) ) / Rb;

    dI = Ia - Ib;

    set_default_plot_settings()
    figure();

    plot(t, za, t, zb)
    xlabel('Time (s)')
    ylabel('z')

    figure();
    plot(t, Ia, t, Ib)
    xlabel('Time (s)')
    ylabel('Current (A)')

    keyboard

end