function out = run_discrete_time_simulation_complete(I_chg, I_dch, ...
               cv_cutoff_current_amps, Qa, Qb, Ra, Rb, za0, zb0, f_ocv, ...
               Vmin, Vmax)
    % Discrete-time simulation of the current imbalance system
    %
    % Simulates a complete charge-discharge cycle, including CC and CV
    % charging operation. Assumes transition from CC to CV charging mode
    % once Vmax is reached.
    %
    % Input charging and discharging current vectors are kept separate
    %
    % Parameters
    % ---------
    %   I_chg:   input current vector in Amperes for charge (+ve is discharge)
    %   I_dch:   input current vector in Amperes for discharge (+ve is discharge)
    %   cv_cutoff_current_amps: cutoff condition for the CV hold (Amperes)
    %   Qa:      capacity of battery A in Ampere-seconds
    %   Qb:      capacity of battery B in Ampere-seconds
    %   Ra:      resistance of battery A in Ohms
    %   Rb:      resistance of battery B in Ohms
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

    R = Ra + Rb;

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
    Ia(1) = (f_ocv(za0) - f_ocv(zb0) + Rb*I_chg(1))/R;
    Ib(1) = (f_ocv(zb0) - f_ocv(za0) + Ra*I_chg(1))/R;
    Vt(1) = f_ocv(za0) - Ia(1)*Ra;
    tt(1) = 0;
    state(1) = 0; % 0: CC charge, 1: CV charge, 2: CC discharge

    % Main loop for CCCV charging
    k = 1;

    t_chg_cc = 0;

    while true

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

        Ia(k+1) = (f_ocv(za(k)) - f_ocv(zb(k)) + Rb*I_chg(k)) / R;
        Ib(k+1) = (f_ocv(zb(k)) - f_ocv(za(k)) + Ra*I_chg(k)) / R;

        Vt(k+1) = f_ocv(za(k)) - Ia(k)*Ra;

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
            I_chg(k+1) = ( f_ocv(za(k))*Rb + f_ocv(zb(k))*Ra - Vmax*(Ra+Rb) ) / (Ra*Rb);

        end

        tt(k+1) = dt*k;
        k = k + 1;

    end



    % Main loop for discharging
    k2 = k;

    while true

        if numel(I_dch) < k2
            error('Not enough data points provided for the current vector')
        end

        if verbose
           fprintf('CC Discharge Step %g: %.3fV\n', k2, Vt(k2))
        end
        
        state(k2+1) = 2;

        Ia(k2+1) = (f_ocv(za(k2)) - f_ocv(zb(k2)) + Rb*I_dch(k2 - k + 1)) / R;
        Ib(k2+1) = (f_ocv(zb(k2)) - f_ocv(za(k2)) + Ra*I_dch(k2 - k + 1)) / R;

        Vt(k2+1) = f_ocv(za(k2)) - Ia(k2)*Ra;

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
