function out = solve_z_dynamics_cccv_complete(tcc, ...
            I_chg, I_dch, I_current_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ...
            ocv_fn, Vmin, Vmax)
    % State of charge dynamics computation including CCCV and subsequent
    % discharge
    %
    % x = [za, zb]^T
    %
    % Parameters
    % ---------
    %   tccv :  time vector (n elements)
    %   I_chg:  current input vector for charge (n elements)
    %   I_dch:  current input vector for discharge (n elements)
    %   I_current_cutoff: CV hold cutoff condition in Amperes
    %   alpha:  OCV by SOC slope parameter, assuming affine OCV function
    %   Ra:     series resistance of cell A
    %   Rb:     series resistance of cell B
    %   Qa:     capacity of cell A in Ampere-seconds
    %   Qb:     capacity of cell B in Ampere-seconds
    %   za0:    initial state of charge for cell A (0,1)
    %   zb0:    initial state of charge for cell B (0,1)
    %   ocv_fn: OCV function = f(z)
    %   Vmin:   target discharge min voltage
    %   Vmax:   target charge max voltage
    %
    % Outputs:
    % ---------
    %   out:    output struct including time vector, Za, Zb, Ia, Ib, and Vt
    %
    % Remarks:
    % I'm not sure how to use lsim to terminate the sim after a certain
    % condition, so I'm going to simulate for a large time vector then 
    % truncate the output.

    %% PART 1: CC Discharge Simulation
    [zacc1, zbcc1, Iacc1, Ibcc1] = solve_cc(tcc, I_chg, Ra, Rb, Qa, Qb, alpha, ...
                                        za0, zb0);

    Vtacc1 = ocv_fn(zacc1) - Iacc1*Ra;
    
    time_at_vmax = interp1(Vtacc1, tcc, Vmax);

    assert(~isnan(time_at_vmax), ...
        'Not enough simulation time given to reach max voltage!')
    
    idx = find(tcc < time_at_vmax);
    zacc1 = zacc1(idx);
    zbcc1 = zbcc1(idx);
    Iacc1 = Iacc1(idx);
    Ibcc1 = Ibcc1(idx);
    Vtacc1 = Vtacc1(idx);
    tcc1 = tcc(idx);

    %% PART II: CV Charge SIMULATION
    tcv = [0:0.1:10*3600]';

    [zacv, zbcv, Iacv, Ibcv] = solve_cv(tcv, Ra, Rb, Qa, Qb, alpha, ...
                                          zacc1(end), zbcc1(end), Vmin, Vmax); 

    idx = find(abs(Iacv + Ibcv) < abs(I_current_cutoff));

    assert(~isempty(idx), ...
          'Not enough simulation time given to reach CV hold cut-off condition.')

    zacv(idx) = [];
    zbcv(idx) = [];
    Iacv(idx) = [];
    Ibcv(idx) = [];
    Vtacv = ones(size(zacv)) * Vmax;
    tcv(idx) = [];

    %% Part III: CC Discharge Simulation
    tcc2 = tcc;
    
    [zacc2, zbcc2, Iacc2, Ibcc2] = solve_cc(tcc2, I_dch, Ra, Rb, Qa, Qb, alpha, ...
                                        zacv(end), zbcv(end));

    Vtacc2 = ocv_fn(zacc2) - Iacc2*Ra;
    
    time_at_vmin = interp1(Vtacc2, tcc, Vmin);

    idx = find(tcc2 < time_at_vmin);
    zacc2 = zacc2(idx);
    zbcc2 = zbcc2(idx);
    Iacc2 = Iacc2(idx);
    Ibcc2 = Ibcc2(idx);
    tcc2 = tcc2(idx);
    Vtacc2 = Vtacc2(idx);

    %% Concatenate the final output vector, combining the CC and CV results
    za = [zacc1 ; zacv ; zacc2];
    zb = [zbcc1 ; zbcv ; zbcc2];
    tfinal = [tcc1 ; tcc1(end) + tcv ; tcc1(end) + tcv(end) + tcc2];
    Ia = [Iacc1 ; Iacv ; Iacc2];
    Ib = [Ibcc1 ; Ibcv ; Ibcc2];
    Vta = [Vtacc1 ; Vtacv ; Vtacc2];

    out.t = tfinal;
    out.za = za;
    out.zb = zb;
    out.Ia = Ia;
    out.Ib = Ib;
    out.Vt = Vta;
    out.t_chg_cc = tcc1(end);
    out.t_chg_cv = tcc1(end) + tcv(end);
    out.t_dch_cc = tfinal(end);

end

function [za, zb, Ia, Ib] = solve_cc(t, I, Ra, Rb, Qa, Qb, alpha, za0, zb0)
    % Solver for the constant current case

    z0 = [za0 zb0]';

    [A, B, C, D] = build_state_matrices_cc_mode(Ra, Rb, Qa, Qb, alpha);
    
    sys = ss(A, B, C, D);
    z = lsim(sys, I, t, z0);

    za = z(:, 1);
    zb = z(:, 2);
 
    % Get the branch currents
    dz = za - zb;

    Ia = (+alpha*dz + I*Rb) / (Ra + Rb);
    Ib = (-alpha*dz + I*Ra) / (Ra + Rb);   

end

function [A, B, C, D] = build_state_matrices_cc_mode(Ra, Rb, Qa, Qb, alpha)

    Alpha = [-alpha, alpha; 
             alpha, -alpha];

    q = [1/Qa 1/Qb]';
    r = [Rb Ra]';
    R = [Ra + Rb];
    A = q ./ R .* Alpha;
    B = - q .* r ./ R;
    C = eye(2);
    D = 0;

end

function [za, zb, Ia, Ib] = solve_cv(t, Ra, Rb, Qa, Qb, alpha, za0, zb0, Vmin, Vmax)
    % Solver for the constant voltage case

    dU = Vmax - Vmin;

    u = ones(size(t)) * dU;
    
    % Solution using lsim
    [A, B, C, D] = build_state_matrices_cv_mode(Ra, Rb, Qa, Qb, alpha);
    sys = ss(A, B, C, D);
    out = lsim(sys, u, t, [za0 ; zb0]);
    
    za = out(:, 1);
    zb = out(:, 2);
    
    taua = Qa*Ra/alpha;
    taub = Qb*Rb/alpha;

    % Analytical form
    za = za0 * exp(-t./taua) - (1/alpha) * (exp(-t./taua) - 1) .* dU;
    zb = zb0 * exp(-t./taub) - (1/alpha) * (exp(-t./taub) - 1) .* dU;

    Ia = (alpha*za - dU) / Ra;
    Ib = (alpha*zb - dU) / Rb;

    
end

function [A, B, C, D] = build_state_matrices_cv_mode(Ra, Rb, Qa, Qb, alpha)

    A = [-alpha/(Qa*Ra)  0               ; ...
         0               -alpha/(Qb*Rb)] ;

    B = [1/(Qa*Ra) ; ...
         1/(Qb*Rb)];

    C = eye(2);
    
    D = 0;

end