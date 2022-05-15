function dI = solve_di_dynamics(t, I, Ra, Rb, Qa, Qb, alpha, dz0)
    % Delta current dynamics
    %
    % dI = Ia - Ib

    R = Ra + Rb;
    kappa = 1/alpha * (Ra*Qa - Rb*Qb) / (Qa + Qb);
    tau = R / alpha * (Qa*Qb) / (Qa + Qb);

    dI = 2 * alpha * dz0 / R * exp(-t./tau) ...
        - (   (2 * alpha * kappa) / R * (exp(-t./tau) - 1) ...
            + (Ra - Rb) / R  ...
          ) .* I;

    % Alternate form
    % dI = alpha / (Ra+Rb) .* dz + - (Ra - Rb) / (Ra + Rb) .* I;

end