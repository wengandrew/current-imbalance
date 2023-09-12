function [condition, zbound_l2, zbound_linf, ibound_l2, ibound_linf] = ...
            solve_imbalance_bounds(ocv, p, t, I, dz0)
    % Compute imbalance bounds based on Hamid's proofs
    % 
    % Args
    % - ocv : an OCV vs z function)
    % - p   : cell parameters struct
    % - t   : time vector with n elements
    % - I   : a vector of currents in Amperes (negative is charge), with n
    %         elements
    % - dz0 : initial soc imbalance = za - zb

    % Outputs:
    % - 
    Ra = p.Ra;
    Rb = p.Rb;
    Qa = p.Qa;
    Qb = p.Qb;

    zz = linspace(0, 1, 1000);
    dOCVdZ = gradient(ocv(zz))./gradient(zz);

    k1 = min(dOCVdZ);
    k2 = max(dOCVdZ);

    A = - (1 / (Ra+Rb)) * (1/Qa + 1/Qb);
    B = abs( 1/(Ra+Rb) * (Rb/Qa - Ra/Qb) );

    signorm_I = sqrt(trapz(t, I.^2));

    condition   = - A * k1 / B;    
    zbound_l2   = - B / (A * k1) * signorm_I + abs(dz0) * sqrt(-1/(2*A*k1));
    zbound_linf = abs(dz0) * exp(k1*A*t) + abs(B/(A*k1))*max(abs(I)) * ( 1 - exp(k1*A*t) ) ;
    
    ibound_l2   = abs( (2*k2*zbound_l2   + (Ra - Rb)*condition) / (Ra+Rb) );
    ibound_linf = abs( (2*k2*zbound_linf + (Ra - Rb)*condition) / (Ra+Rb) );
    

end

%     signorm_z_tilde = sqrt(sum(z_tilde.^2));
%     signorm_I = sqrt(sum(I.^2));

%     signorm_I = norm(I, 2);
%     signorm_z_tilde = norm(z_tilde, 2);
