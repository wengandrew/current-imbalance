function p = initialize_cell_props_dsmc()
    % Set default cell properties here

    p.Qa = 3 * 3600; % Amp-seconds
    p.Ra = 150 / 1000; % Ohms

    q = 0.7;
    r = 1.1;

    p.Qb = p.Qa/q;
    p.Rb = p.Ra/r;

    p.U0 = 3.0;
    p.alpha = 1.2;

    p.f_Ua = @(z) U0 + alpha * z;
    p.f_Ub = @(z) U0 + alpha * z;

end
