function p = initialize_cell_props()
    % Set default cell properties here

    p.Qa = 4 * 3600; % Amp-seconds
    p.Qb = 5 * 3600; % Amp-seconds
    p.Ra = 60 / 1000; % Ohms
    p.Rb = 50 / 1000; % Ohms

    p.U0 = 3.0;
    p.alpha = 1.2;

    p.f_Ua = @(z) U0 + alpha * z;
    p.f_Ub = @(z) U0 + alpha * z;
    
end