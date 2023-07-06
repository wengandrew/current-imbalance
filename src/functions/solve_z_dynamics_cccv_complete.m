function out = solve_z_dynamics_cccv_complete(tcc, ...
            I_chg, I_dch, I_current_cutoff, Ra, Rb, Qa, Qb, za0, zb0, ...
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
    %   Ra:     series resistance of cell A
    %   Rb:     series resistance of cell B
    %   Qa:     capacity of cell A in Ampere-seconds
    %   Qb:     capacity of cell B in Ampere-seconds
    %   za0:    initial state of charge for cell A (0,1)
    %   zb0:    initial state of charge for cell B (0,1)
    %   ocv_fn: OCV function = f(z), z in (0,1)
    %
    % Outputs:
    % ---------
    %   out:    output struct including time vector, Za, Zb, Ia, Ib, and Vt
    %
    % Remarks:
    % I'm not sure how to use lsim to terminate the sim after a certain
    % condition, so I'm going to simulate for a large time vector then 
    % truncate the output.

    alpha = ocv_fn(1) - ocv_fn(0);

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
    out.za_ohmic = [cc1.za_ohmic ; cv1.za_ohmic ; cc2.za_ohmic];
    out.za_rebal = [cc1.za_rebal ; cv1.za_rebal ; cc2.za_rebal];
    out.zb_ohmic = [cc1.zb_ohmic ; cv1.zb_ohmic ; cc2.zb_ohmic];
    out.zb_rebal = [cc1.zb_rebal ; cv1.zb_rebal ; cc2.zb_rebal];
    out.Ia = [cc1.Ia ; cv1.Ia ; cc2.Ia];
    out.Ib = [cc1.Ib ; cv1.Ib ; cc2.Ib];
    out.Ia_ohmic = [cc1.Ia_ohmic ; cv1.Ia_ohmic ; cc2.Ia_ohmic];
    out.Ia_rebal = [cc1.Ia_rebal ; cv1.Ia_rebal ; cc2.Ia_rebal];
    out.Ib_ohmic = [cc1.Ib_ohmic ; cv1.Ib_ohmic ; cc2.Ib_ohmic];
    out.Ib_rebal = [cc1.Ib_rebal ; cv1.Ib_rebal ; cc2.Ib_rebal];
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

            za_rebal = 1/Q * ( (Qa + Qb*exp(-t/tau))*za0 + Qb*(1-exp(-t/tau))*zb0 );
            za_ohmic = 1/Q * ( +Qb*kappa*(1-exp(-t/tau)) - t) * u;
            za = za_rebal + za_ohmic;

            zb_rebal = 1/Q * ( (Qa*(1-exp(-t/tau)))*za0 + (Qb + Qa*exp(-t/tau))*zb0 );
            zb_ohmic = 1/Q * ( -Qa*kappa*(1-exp(-t/tau)) - t) * u;
            zb = zb_rebal + zb_ohmic;
            
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
    out.za_rebal = za_rebal;
    out.za_ohmic = za_ohmic;
    out.zb_rebal = zb_rebal;
    out.zb_ohmic = zb_ohmic;
    out.Ia = Ia;
    out.Ib = Ib;
    out.Ia_ohmic = Ia_ohmic;
    out.Ia_rebal = Ia_rebal;
    out.Ib_ohmic = Ib_ohmic;
    out.Ib_rebal = Ib_rebal;

    out = struct2table(out);

    % Truncate data past Vlim
    time_at_vlim = interp1(Vt(~isnan(Vt)), t(~isnan(Vt)), Vlim);
    assert(~isnan(time_at_vlim), ...
        'Not enough simulation time given to reach voltage limit.')

    idx = find(t >= time_at_vlim);
    out(idx, :) = [];

end


function out = solve_cv(t, R2, R1, Q2, Q1, alpha, z20, z10, Vmin, Vmax, I_current_cutoff)
    % Solver for the constant voltage case

    Z_SOLUTION_METHOD = 'analytic';

    dU = Vmax - Vmin;
    u = ones(size(t)) * dU;

    Q = Q1 + Q2;
    kappa = (R2*Q2-R1*Q1)/(Q1+Q2)/alpha;
    tau = (R1+R2)/alpha*(Q1*Q2)/(Q1+Q2);

    switch Z_SOLUTION_METHOD

        case 'lsim'

            [A, B, C, D] = build_state_matrices_cv_mode(R2, R1, Q2, Q1, alpha);
            sys = ss(A, B, C, D);
            sim_out = lsim(sys, u, t, [z20 ; z10]);
            
            za = sim_out(:, 1);
            zb = sim_out(:, 2);

        case 'analytic'

            taua = Q2*R2/alpha;
            taub = Q1*R1/alpha;
        
            za = z20 * exp(-t./taua) - (1/alpha) * (exp(-t./taua) - 1) .* dU;
            zb = z10 * exp(-t./taub) - (1/alpha) * (exp(-t./taub) - 1) .* dU;

    end

    Ia1 = alpha*za / R2;
    Ia2 = -u / R2;

    Ib1 = alpha*zb / R1;
    Ib2 = -u / R1;

    Ia = Ia1 + Ia2;
    Ib = Ib1 + Ib2;

    Vt = ones(size(za)) * Vmax;
    
    out.t = t;
    out.Vt = Vt;
    out.za = za;
    out.zb = zb;
    out.za_ohmic = 1/Q * ( +Q1*kappa*(1-exp(-t./tau)) - t) .* (Ia + Ib);
    out.za_rebal = 1/Q * ( (Q2 + Q1*exp(-t./tau))*z20 + Q1*(1-exp(-t./tau))*z10 );
    out.zb_ohmic = 1/Q * ( -Q2*kappa*(1-exp(-t./tau)) - t) .* (Ia + Ib);
    out.zb_rebal = 1/Q * ( (Q2*(1-exp(-t./tau)))*z20 + (Q1 + Q2*exp(-t./tau))*z10 );
    out.Ia = Ia;
    out.Ib = Ib;
    out.Ia_ohmic = R1/(R1+R2) * (Ia+Ib);
    out.Ia_rebal = Ia - out.Ia_ohmic;
    out.Ib_ohmic = R2/(R1+R2) * (Ia+Ib);
    out.Ib_rebal = Ib - out.Ib_ohmic;

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