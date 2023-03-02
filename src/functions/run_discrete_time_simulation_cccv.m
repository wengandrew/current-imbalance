function out = run_discrete_time_simulation_cccv(t, I, Qa, Qb, Ra, Rb, ...
    za0, zb0, ocv, Vmax)
    % Discrete-time simulation of the current imbalance system
    %
    % Supports both CC and CV operation. Assumes transition from CC to CV
    % charging mode once Vmax is reached.
    %
    % Parameters
    % ---------
    %   t: input time vector in seconds
    %   I: input current vector in Amperes (+ve is discharge)
    %   Qa: capacity of battery A in Ampere-seconds
    %   Qb: capacity of battery B in Ampere-seconds
    %   Ra: resistance of battery A in Ohms
    %   Rb: resistance of battery B in Ohms
    %   za0: initial SOC of battery A
    %   zb0: initial SOC of battery B
    %   ocv: OCV function = f(z); an anonymous function
    %   Vmax: voltage termination condition
    %
    % Assumptions
    % ---------
    %  - fixed time-step simulation
    %  - charge simulation ends on Vmax
    %
    % The simulation does not assume constant current

    dt = t(2) - t(1);

    R = Ra + Rb;

    % Initial conditions
    za = zeros(size(t));
    zb = zeros(size(t));
    Ia = zeros(size(t));
    Ib = zeros(size(t));
    Vt = zeros(size(t));

    za(1) = za0;
    zb(1) = zb0;
    Ia(1) = (ocv(za0) - ocv(zb0) + Rb*I(1))/R;
    Ib(1) = (ocv(zb0) - ocv(za0) + Ra*I(1))/R;
    Vt(1) = ocv(za0) - Ia(1)*Ra;
    state(1) = 0; % 0: Constant Current, 1: Constant Voltage

    % Main loop
    for k = 1:numel(t) - 1

        za(k+1) = za(k) - dt/Qa * Ia(k);
        zb(k+1) = zb(k) - dt/Qb * Ib(k);

        za(k+1) = min(za(k+1), 1);
        zb(k+1) = min(zb(k+1), 1);

        Ia(k+1) = (ocv(za(k)) - ocv(zb(k)) + Rb*I(k)) / R;
        Ib(k+1) = (ocv(zb(k)) - ocv(za(k)) + Ra*I(k)) / R;

        Vt(k+1) = ocv(za(k)) - Ia(k)*Ra;

        % Check for potentiostatic mode, in which case update current
        if Vt(k+1) >= Vmax
            Vt(k+1) = Vmax;
            I(k+1) = ( ocv(za(k))*Rb + ocv(zb(k))*Ra - Vmax*(Ra+Rb) ) / (Ra*Rb);
        end

        assert(~isnan(Vt(k+1)), 'Terminal voltage is NaN.')

    end

    assert(~any(isnan(Vt)), 'NaN found in terminal voltage. Check OCV function.')

    % Assemble output as a struct
    out.t = t;
    out.za = za;
    out.zb = zb;
    out.Ia = Ia;
    out.Ib = Ib;
    out.Vt = Vt;

end
