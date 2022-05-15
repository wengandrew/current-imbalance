function dz = solve_dz_dynamics(t, I, alpha, Ra, Rb, Qa, Qb, dz0)
    % Delta z dynamics computation
    %
    % Delta z = za - zb
    %
    % Parameters
    % ---------
    %   t:     time input vector (n elements)
    %   I:     current input vector (n elements)
    %   alpha: OCV by SOC slope parameter, assuming affine OCV function
    %   Ra:    series resistance of cell A
    %   Rb:    series resistance of cell B
    %   Qa:    capacity of cell A in Ampere-seconds
    %   Qb:    capacity of cell B in Ampere-seconds
    %   dz0:   initial za - zb (number between 0 to 1)
    %
    % Outputs
    % --------
    %   dz: output vector (n elements)
    
    % State space formulation
    %     A = -alpha/(Ra+Rb) * (1/Qa + 1/Qb);
    %     B = (Ra/Qb - Rb/Qa) / (Ra+Rb);
    %     C = 1;
    %     D = 0;

    %     dz = dz0 * exp(A.*t) + B/A*(exp(A*t) - 1).*I;
    %
    % Alternatively, use the built-in solver:
    % sys = ss(A, B, C, D);
    % dz = lsim(sys, I, t, dz0);

    % Define A = -tau^-1
    % Define B = kappa * tau^-1

    R = Ra + Rb;

    tau = R / alpha * (Qa*Qb) / (Qa + Qb);
    kappa = 1 / alpha * (Ra*Qa - Rb*Qb) / (Qa + Qb);

    dz = dz0 * exp(-t ./ tau) - kappa * (exp(-t./tau) - 1) .* I;

end
