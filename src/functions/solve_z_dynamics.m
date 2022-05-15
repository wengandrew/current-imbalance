function [za, zb] = solve_z_dynamics(t, I, alpha, Ra, Rb, Qa, Qb, za0, zb0)
    % State of charge dynamics computation
    %
    % x = [za, zb]^T
    %
    % Parameters
    % ---------
    %   t:      time input vector (n elements)
    %   I:      current input vector (n elements)
    %   alpha:  OCV by SOC slope parameter, assuming affine OCV function
    %   Ra:     series resistance of cell A
    %   Rb:     series resistance of cell B
    %   Qa:     capacity of cell A in Ampere-seconds
    %   Qb:     capacity of cell B in Ampere-seconds
    %   za0:    initial state of charge for cell A (0,1)
    %   zb0:    initial state of charge for cell B (0,1)
    %
    % Outputs:
    % ---------
    %   za:     output vector of state of charge for cell A (n elements)
    %   zb:     output vector of state of charge for cell B (n elements)

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
    z = lsim(sys, I, t, z0);

    za = z(:, 1);
    zb = z(:, 2);
 
end