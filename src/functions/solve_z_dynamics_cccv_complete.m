function out = solve_z_dynamics_cccv_complete(tcc, ...
            I_chg, I_dch, I_current_cutoff, alpha, Ra, Rb, Qa, Qb, za0, zb0, ...
            ocv_fn, Vmin, Vmax)
    % State of charge dynamics computation for a full charge-discharge
    % cycle, including:
    %
    % 1. Constant current charge
    % 2. Constant voltage charge
    % 3. Constant current discharge
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
    %   ocv_fn: OCV function = f(z), z in (0,1)
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

    % Constant current charge
    cc1 = solve_cc(tcc, I_chg, Ra, Rb, Qa, Qb, alpha, za0, zb0, ocv_fn, Vmax);

    % Constant voltage charge
    tcv = [0:0.1:10*3600]';
    cv1 = solve_cv(tcv, Ra, Rb, Qa, Qb, alpha, cc1.za(end), cc1.zb(end), Vmin, Vmax, I_current_cutoff); 

    % Constant current discharge
    tcc2 = tcc;
    cc2 = solve_cc(tcc2, I_dch, Ra, Rb, Qa, Qb, alpha, cv1.za(end), cv1.zb(end), ocv_fn, Vmin);

    % Combine results
    out.t = [cc1.t ; cc1.t(end) + cv1.t ; cc1.t(end) + cv1.t(end) + cc2.t];
    out.Vt = [cc1.Vt ; cv1.Vt ; cc2.Vt];
    out.za = [cc1.za ; cv1.za ; cc2.za];
    out.zb = [cc1.zb ; cv1.zb ; cc2.zb];
    out.Ia = [cc1.Ia ; cv1.Ia ; cc2.Ia];
    out.Ib = [cc1.Ib ; cv1.Ib ; cc2.Ib];
    out.t_chg_cc = cc1.t(end);
    out.t_chg_cv = cc1.t(end) + cv1.t(end);
    out.t_dch_cc = out.t(end);

end

function out = solve_cc(t, I, Ra, Rb, Qa, Qb, alpha, za0, zb0, ocv_fn, Vlim)
    %% Solver for the constant current case

    Z_SOLUTION_METHOD = 'analytic';

    R = Ra + Rb;

    z0 = [za0 zb0]';

    switch Z_SOLUTION_METHOD

        case 'lsim'

            [A, B, C, D] = build_state_matrices_cc_mode(Ra, Rb, Qa, Qb, alpha);

            sys = ss(A, B, C, D);
            z = lsim(sys, I, t, z0);
        
            za = z(:, 1);
            zb = z(:, 2);

        case 'analytic'

            Q = Qa + Qb;
            u = I(1);

            tau = (Ra + Rb) / alpha * (Qa * Qb) / Q;
            kappa = 1/alpha * (Ra*Qa - Rb*Qb) / (Qa + Qb);

            za = 1/Q * ( (Qa + Qb*exp(-t/tau))*za0 + Qb*(1-exp(-t/tau))*zb0 ) + ...
                 1/Q * ( +Qb*kappa*(1-exp(-t/tau)) - t) * u;

            zb = 1/Q * ( (Qa*(1-exp(-t/tau)))*za0 + (Qb + Qa*exp(-t/tau))*zb0 ) + ...
                 1/Q * ( -Qa*kappa*(1-exp(-t/tau)) - t) * u;
            
    end
  
    % Get the branch currents
    dz = za - zb;

    Ia_ohmic = I*Rb / R;
    Ia_rebal = +alpha*dz / R;
    Ib_ohmic = I*Ra / R;
    Ib_rebal = -alpha*dz / R;

    Ia = Ia_ohmic + Ia_rebal;
    Ib = Ib_ohmic + Ib_rebal;

    Vt = ocv_fn(za) - Ia*Ra;

    out.t = t;
    out.Vt = Vt;
    out.za = za;
    out.zb = zb;
    out.Ia = Ia;
    out.Ib = Ib;
    out.Ia_ohmic = Ia_ohmic;
    out.Ia_rebal = Ia_rebal;
    out.Ib_ohmic = Ib_ohmic;
    out.Ib_rebal = Ib_rebal;

    out = struct2table(out);

    % Truncate data past Vlim
    time_at_vlim = interp1(Vt, t, Vlim);
    assert(~isnan(time_at_vlim), ...
        'Not enough simulation time given to reach voltage limit.')

    idx = find(t >= time_at_vlim);
    out(idx, :) = [];

end


function out = solve_cv(t, Ra, Rb, Qa, Qb, alpha, za0, zb0, Vmin, Vmax, I_current_cutoff)
    % Solver for the constant voltage case

    Z_SOLUTION_METHOD = 'analytic';

    dU = Vmax - Vmin;
    u = ones(size(t)) * dU;
    
    switch Z_SOLUTION_METHOD

        case 'lsim'

            [A, B, C, D] = build_state_matrices_cv_mode(Ra, Rb, Qa, Qb, alpha);
            sys = ss(A, B, C, D);
            sim_out = lsim(sys, u, t, [za0 ; zb0]);
            
            za = sim_out(:, 1);
            zb = sim_out(:, 2);

        case 'analytic'

            taua = Qa*Ra/alpha;
            taub = Qb*Rb/alpha;
        
            za = za0 * exp(-t./taua) - (1/alpha) * (exp(-t./taua) - 1) .* dU;
            zb = zb0 * exp(-t./taub) - (1/alpha) * (exp(-t./taub) - 1) .* dU;

    end

    Ia1 = alpha*za / Ra;
    Ia2 = -u / Ra;

    Ib1 = alpha*zb / Rb;
    Ib2 = -u / Rb;

    Ia = Ia1 + Ia2;
    Ib = Ib1 + Ib2;

    Vt = ones(size(za)) * Vmax;
    
    out.t = t;
    out.Vt = Vt;
    out.za = za;
    out.zb = zb;
    out.Ia = Ia;
    out.Ib = Ib;
    out.Ia1 = Ia1;
    out.Ia2 = Ia2;
    out.Ib1 = Ib1;
    out.Ib2 = Ib2;

    out = struct2table(out);

    % Truncate data past cutoff-current
    idx = find(abs(Ia + Ib) <= abs(I_current_cutoff));
    assert(~isempty(idx), ...
          'Not enough simulation time provided to reach CV hold cut-off condition.')
    out(idx, :) = [];

end


function [A, B, C, D] = build_state_matrices_cc_mode(Ra, Rb, Qa, Qb, alpha)

    Alpha = [-alpha, +alpha; 
             +alpha, -alpha];

    q = [1/Qa 1/Qb]';
    r = [Rb Ra]';
    R = [Ra + Rb];
    A = q ./ R .* Alpha;
    B = - q .* r ./ R;
    C = eye(2);
    D = 0;

end


function [A, B, C, D] = build_state_matrices_cv_mode(Ra, Rb, Qa, Qb, alpha)

    A = [-alpha/(Qa*Ra)  0               ; ...
         0               -alpha/(Qb*Rb)] ;

    B = [1/(Qa*Ra) ; ...
         1/(Qb*Rb)];

    C = eye(2);
    
    D = 0;

end