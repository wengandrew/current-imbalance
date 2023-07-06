function out = run_discrete_time_simulation_complete_nonlinr(I_chg, I_dch, ...
               cv_cutoff_current_amps, Qa, Qb, Ra, Rb, za0, zb0, f_ocva, f_ocvb, ...
               Vmin, Vmax)
    % Discrete-time simulation of the current imbalance system
    %
    % Simulates a complete charge-discharge cycle, including CC and CV
    % charging operation. Assumes transition from CC to CV charging mode
    % once Vmax is reached.
    %
    % Input charging and discharging current vectors are kept separate
    %
    % "nonlinr" version 7/3/2023 includes support for:
    % - two nonlinear OCV curves, one for each cell
    % - Ra and Rb become functions of SOC
    % - implemented a rest step at the end of CCCV charge to match experimental data
    %
    % This function is used mainly to run the model vs experiment comparison.
    % We are keeping this code separate for now since this code supports more
    % functionality like inclusion of rest steps, multiple non-linear OCV
    % curves, and SOC-dependent resistances. This runner can be made part of
    % the main code in a future revision.
    % 
    % Parameters
    % ---------
    %   I_chg:   input current vector in Amperes for charge (+ve is discharge)
    %   I_dch:   input current vector in Amperes for discharge (+ve is discharge)
    %   cv_cutoff_current_amps: cutoff condition for the CV hold (Amperes)
    %   Qa:      capacity of battery A in Ampere-seconds
    %   Qb:      capacity of battery B in Ampere-seconds
    %   Ra:      resistance function for battery A in Ohms
    %   Rb:      resistance function for battery B in Ohms
    %   za0:     initial SOC of battery A
    %   zb0:     initial SOC of battery B
    %   f_ocv:   OCV function = f(z); an anonymous function
    %   Vmax:    voltage termination condition on charge
    %   Vmin:    voltage termination condition on discharge
    %
    % Outputs
    % --------
    %   out:     a struct holding simulation outputs
    %
    %
    % Assumptions
    % ---------
    %  - fixed time-step simulation
    %  - charge simulation ends on Vmax
    %
    % The simulation does not need to assume constant current input

    dt = 1;

    verbose = false;

    % Initialize arrays and initial conditions. We don't how big the arrays
    % will need to be so we just allocate a little bit to start with. At
    % least we enforce the right dimensions.
    za = zeros(10, 1);
    zb = zeros(10, 1);
    Ia = zeros(10, 1);
    Ib = zeros(10, 1);
    Vt = zeros(10, 1);
    tt = zeros(10, 1);
    state = zeros(10, 1);

    za(1) = za0;
    zb(1) = zb0;
    Ia(1) = (f_ocva(za0) - f_ocvb(zb0) + Rb(zb0)*I_chg(1))/(Ra(za0) + Rb(zb0));
    Ib(1) = (f_ocvb(zb0) - f_ocva(za0) + Ra(za0)*I_chg(1))/(Ra(za0) + Rb(zb0));
    Vt(1) = f_ocva(za0) - Ia(1)*Ra(za0);
    tt(1) = 0;
    state(1) = 0; % 0: CC charge, 1: CV charge, 2: CC discharge

    % Main loop for CCCV charging
    k = 1;

    t_chg_cc = 0;

    while true

        tt(k+1) = dt*k;

        if verbose
            fprintf('CC Charge Step %g: %.3f V\n', k, Vt(k))
        end

        state(k+1) = 0;
        za(k+1) = za(k) - dt/Qa * Ia(k);
        zb(k+1) = zb(k) - dt/Qb * Ib(k);

        if za(k+1) > 1.5
            keyboard
        end
        if zb(k+1) > 1.5
            keyboard
        end

        Ia(k+1) = (f_ocva(za(k)) - f_ocvb(zb(k)) + Rb(zb(k))*I_chg(k)) / (Ra(za(k)) + Rb(zb(k)));
        Ib(k+1) = (f_ocvb(zb(k)) - f_ocva(za(k)) + Ra(za(k))*I_chg(k)) / (Ra(za(k)) + Rb(zb(k)));

        Vt(k+1) = f_ocva(za(k)) - Ia(k)*Ra(za(k));

        % Check for potentiostatic (CV) mode, in which case update current
        % while holding potential constant

        if Vt(k+1) >= Vmax

            % Update time counters
            if t_chg_cc == 0
                t_chg_cc = dt*k;
            end
            t_chg_cv = dt*k;

            if verbose
                fprintf('CV Charge Step %g: %.3f V\n', k, Vt(k))
            end

            % CV charge termination condition
            if abs(I_chg(k)) < abs(cv_cutoff_current_amps)
                break
            end

            state(k+1) = 1;
            I_chg(k+1) = ( f_ocva(za(k))*Rb(zb(k)) + f_ocvb(zb(k))*Ra(za(k)) - ...
                Vmax*(Ra(za(k))+Rb(zb(k))) ) / (Ra(za(k))*Rb(zb(k)));

        end

        k = k + 1;

    end

    % Main loop for resting
    kr = k;
    
    REST_DURATION_HRS = 0.5;

    while (kr - k) * dt < REST_DURATION_HRS * 3600 % Rest time

        if verbose
           fprintf('Rest Step %g: %.3fV\n', kr, Vt(kr))
        end
        
        state(kr+1) = 3;

        Ia(kr+1) = (f_ocva(za(kr)) - f_ocvb(zb(kr)) ) / (Ra(za(kr)) + Rb(zb(kr)));
        Ib(kr+1) = (f_ocvb(zb(kr)) - f_ocva(za(kr)) ) / (Ra(za(kr)) + Rb(zb(kr)));

        Vt(kr+1) = f_ocva(za(kr));

        za(kr+1) = za(kr);
        zb(kr+1) = zb(kr);

        tt(kr+1) = dt*kr;

        kr = kr + 1;

    end

    % Main loop for discharging
    k2 = kr;

    while true

        if numel(I_dch) < k2
            error('Not enough data points provided for the current vector')
        end

        if verbose
           fprintf('CC Discharge Step %g: %.3fV\n', k2, Vt(k2))
        end
        
        state(k2+1) = 2;

        Ia(k2+1) = (f_ocva(za(k2)) - f_ocvb(zb(k2)) + Rb(zb(k2))*I_dch(k2 - k + 1)) / (Ra(za(k2)) + Rb(zb(k2)));
        Ib(k2+1) = (f_ocvb(zb(k2)) - f_ocva(za(k2)) + Ra(za(k2))*I_dch(k2 - k + 1)) / (Ra(za(k2)) + Rb(zb(k2)));

        Vt(k2+1) = f_ocva(za(k2)) - Ia(k2)*Ra(za(k2));

        za(k2+1) = za(k2) - dt/Qa * Ia(k2);
        zb(k2+1) = zb(k2) - dt/Qb * Ib(k2);

        tt(k2+1) = dt*k2;

        % Terminate discharge when voltage reaches Vmin
        if Vt(k2+1) <= Vmin
            break
        end

        k2 = k2 + 1;

    end

    assert(~any(isnan(Vt)), ...
        'NaN found in the terminal voltage. Check the OCV function.')

    % Assemble output as a struct
    out.state = state;
    out.t = tt;
    out.t_chg_cc = t_chg_cc;
    out.t_chg_cv = t_chg_cv;
    out.za = za;
    out.zb = zb;
    out.Ia = Ia;
    out.Ib = Ib;
    out.Vt = Vt;

end
