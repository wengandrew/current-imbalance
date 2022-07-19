function [tfinal, za, zb, Ia, Ib] = solve_z_dynamics_cccv(tcv, I, alpha, Ra, Rb, Qa, Qb, za0, zb0, Uf, U0)
    % State of charge dynamics computation including CCCV
    %
    % x = [za, zb]^T
    %
    % Parameters
    % ---------
    %   t:      time input vector (m elements)
    %   I:      current input vector (n elements)
    %   alpha:  OCV by SOC slope parameter, assuming affine OCV function
    %   Ra:     series resistance of cell A
    %   Rb:     series resistance of cell B
    %   Qa:     capacity of cell A in Ampere-seconds
    %   Qb:     capacity of cell B in Ampere-seconds
    %   za0:    initial state of charge for cell A (0,1)
    %   zb0:    initial state of charge for cell B (0,1)
    %   Uf:     target final max voltage
    %   U0:     initial voltage (at SOC = 0)
    %
    % Outputs:
    % ---------
    %   tfinal: output time
    %   za:     output vector of state of charge for cell A (n elements)
    %   zb:     output vector of state of charge for cell B (n elements)
    %   Ia:     output vector of current for cell A (n elements)
    %   Ib:     output vector of current for cell B (n elements)

    %% Simulate the CC portion
    Alpha = [-alpha, alpha; 
             alpha, -alpha];

    q = [1/Qa 1/Qb]';
    r = [Rb Ra]';
    R = [Ra + Rb];
    A = q ./ R .* Alpha;
    B = - q .* r ./ R;
    C = eye(2);
    D = 0;

    z0 = [za0 zb0]';

    sys = ss(A, B, C, D);
    z = lsim(sys, I, tcv, z0);

    zacc = z(:, 1);
    zbcc = z(:, 2);
 
    % Get the branch currents
    dzcc = zacc - zbcc;

    Iacc = (+alpha*dzcc + I*Rb) / (Ra + Rb);
    Ibcc = (-alpha*dzcc + I*Ra) / (Ra + Rb);

    Vta = U0 + alpha.*zacc - Iacc*Ra;
    time_at_vmax = interp1(Vta, tcv, Uf);

    if isempty(time_at_vmax)
        error 'Not enough simulation time given to reach max voltage!'
    end
    
    idx = find(tcv < time_at_vmax);
    zacc = zacc(idx);
    zbcc = zbcc(idx);
    Iacc = Iacc(idx);
    Ibcc = Ibcc(idx);
    tcc = tcv(idx);

    %% Simulate the CV portion
    za0 = zacc(end);
    zb0 = zbcc(end);

    A = [-alpha/(Qa*Ra) 0 ; ...
         0 -alpha/(Qb*Rb)];

    B = [1/(Qa*Ra) ; ...
         1/(Qb*Rb)];
    C = eye(2);
    D = 0;

    tcv = [0:0.1:5*3600]';
    u = ones(size(tcv)) * (Uf - U0);

    sys = ss(A, B, C, D);
    out = lsim(sys, u, tcv, [za0 ; zb0]);

    zacv = out(:, 1);
    zbcv = out(:, 2);

    Iacv = (alpha*zacv - (Uf - U0)) / Ra;
    Ibcv = (alpha*zbcv - (Uf - U0)) / Rb;

    %% Concatenate the final output vector, combining the CC and CV results
    % together
    
    za = [zacc ; zacv];
    zb = [zbcc ; zbcv];
    tfinal = [tcc ; tcv + tcc(end)];
    Ia = [Iacc ; Iacv];
    Ib = [Ibcc ; Ibcv];

end