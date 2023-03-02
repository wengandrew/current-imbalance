function aw20221221_imbalance_bound_simulations(analysis_type, q, r, current_target)
% 
% 
% Args:
%  analysis_type : 'lfp', 'nmc', 'quad'ratic function
%  q : Qa/Qb
%  r : Ra/Rb
%  current_target : target current in amperes

    to_plot = true;

    set_default_plot_settings_manuscript()

    switch analysis_type
        case 'lfp'
            alpha = 0.6;
            Vmax = 3.6;
            ocv = load_ocv_fn('lfp');
        case 'nmc'
            alpha = 1.2;
            Vmax = 4.2;
            ocv = load_ocv_fn('nmc');
        case 'quad'
            alpha = 1.2;
            Vmax = 4.2;
            ocv = @(x) 0.6.*x.^2 + 0.6.*x + 3;
    end

    % Initialize model parameters
    Ra = 0.05; % ohms
    Qa = 5 * 3600; % As
    alpha_nmc = 1.2;
    alpha_lfp = 0.6;
    za0 = 0.00;
    zb0 = 0.00;
    U0  = 3.0;
    Vmin = 3.0;

    Qb = Qa/q;
    Rb = Ra/r;
        
    zz = linspace(0, 1, 1000);
%     current_target = -Qa / (3 * 3600);

    %% Check current imbalance conditions
    dOCVdZ = gradient(ocv(zz))./gradient(zz);

    if to_plot

        % Inspect the OCV function
        figure()
        ax1 = subplot(211);
        plot(zz, ocv(zz));
        xlabel('z')
        ylabel('OCV')
        xlim([0, 1])
    
        ax2 = subplot(212);
        plot(zz, dOCVdZ);
        xlabel('z')
        ylabel('dOCV/dz')
        xlim(ax1, [0, 1])
        linkaxes([ax1, ax2], 'x')

    end

    %% Test the analytic solution
    % Initialize simulation vectors
    t_max_hours = 6 * 2.5/abs(current_target);
    t = linspace(0, t_max_hours*3600, 5e4)';
%     t = [0:60:t_max_hours*3600]';
%     keyboard
    I = current_target*ones(size(t)); % applied current (A)
%     I = current_target .* linspace(1, 0, numel(t))' .* 1 .* sin(t);
    [tfinal, za, zb, Ia, Ib] = solve_z_dynamics_cccv(t, I, ...
        alpha, Ra, Rb, Qa, Qb, za0, zb0, Vmax, U0);

    %% Test the simulation based solution
    out = run_discrete_time_simulation_cccv(t, I, Qa, Qb, Ra, Rb, ...
            za0, zb0, ocv, Vmax);

    %% Current imbalance bound evaluation from Hamid
    k1 = min(dOCVdZ);
    k2 = max(dOCVdZ);
    A = - (1 / (Ra+Rb)) * (1/Qa + 1/Qb);
    B = abs( 1/(Ra+Rb) * (Rb/Qa - Ra/Qb) );

    kappa = 1 / alpha * (Ra*Qa - Rb*Qb) / (Qa+Qb);

    fprintf('Qa = %.1f Ah\n', Qa/3600)
    fprintf('Ra = %.1f mOhms\n', Ra*1000)
    fprintf('Qa/Qb = %.1f\n', Qa/Qb )
    fprintf('Ra/Rb = %.1f\n', Ra/Rb )
    fprintf('I = %.1f A \n', max(abs(I)))

    z_tilde = out.zb - out.za;
    z0_tilde = zb0 - za0;

    signorm_z_tilde = sqrt(trapz(out.t, z_tilde.^2));
    signorm_I = sqrt(trapz(out.t, I.^2));
% 
%     signorm_z_tilde = sqrt(sum(z_tilde.^2));
%     signorm_I = sqrt(sum(I.^2));

%     signorm_I = norm(I, 2);
%     signorm_z_tilde = norm(z_tilde, 2);

    L2condition_LHS = abs(current_target);
    L2condition_RHS = - A * k1 / B;
    
    L2bound_LHS = signorm_z_tilde;
    L2bound_RHS = - B / (A * k1) * signorm_I + abs(z0_tilde) * sqrt(-1/(2*A*k1));

    Linfcondition_LHS = max(abs(I));
    Linfcondition_RHS = - A * k1 / B;

    Linfbound_LHS = max(z_tilde);
    Linfbound_RHS = - B / (A * k1) * max(abs(I)) + abs(z0_tilde);

    fprintf('L2   Condition:  %.3f <= %.3f ? \t%g \n', L2condition_LHS, L2condition_RHS, L2condition_LHS < L2condition_RHS);
    fprintf('L2   Bound:      %.3f <= %.3f ? \t%g \n', L2bound_LHS, L2bound_RHS, L2bound_LHS < L2bound_RHS)
    
    fprintf('Linf Condition:  %.3f <= %.3f ? \t%g \n', Linfcondition_LHS, Linfcondition_RHS, Linfcondition_LHS < Linfcondition_RHS);
    fprintf('Linf Bound:      %.3f <= %.3f ? \t%g \n', Linfbound_LHS, Linfbound_RHS, Linfbound_LHS < Linfbound_RHS);

    % See the results

    if to_plot

        figure('OuterPosition', [0, 0, 900, 1200])

        ax3 = subplot(411);
        line(tfinal./3600, za, 'Color', 'b', 'DisplayName', 'Analytic')
        line(tfinal./3600, zb, 'Color', 'r', 'DisplayName', 'Analytic')
        line(out.t./3600, out.za, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
        line(out.t./3600, out.zb, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
        xlabel('Time (hrs)')
        ylabel('$z$', 'Interpreter', 'Latex')
        legend show

        ax4 = subplot(412);
        line(tfinal./3600, abs(za - zb), 'Color', 'r', 'DisplayName', 'Analytic')
        yline(abs(kappa * current_target), 'DisplayName', 'Bound ($\kappa I$)', 'LineStyle', '--', 'Color', 'r')
        line(out.t./3600, abs(out.za - out.zb), 'Color', 'b', 'DisplayName', 'Simulated')
        yline(Linfbound_RHS, 'DisplayName', 'Bound ($-\frac{B}{Ak_1}\max(|I|) + |\tilde{z}_0|$)', 'LineStyle', '--', 'Color', 'b')

        xlabel('Time (hrs)')
        ylabel('$|\Delta z|$', 'Interpreter', 'Latex')
        legend show
    
        ax5 = subplot(413);
        line(tfinal./3600, Ia, 'Color', 'b', 'DisplayName', 'Analytic')
        line(tfinal./3600, Ib, 'Color', 'r', 'DisplayName', 'Analytic')
        line(out.t./3600, out.Ia, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
        line(out.t./3600, out.Ib, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
        xlabel('Time (hrs)')
        ylabel('Current (A)')
        legend show
    
        ax6 = subplot(414);
        line(tfinal/3600, U0 + alpha*za - Ia*Ra, 'Color', 'k', 'DisplayName', 'Analytic')
        line(out.t./3600, out.Vt, 'Color', 'k', 'LineStyle', ':', 'DisplayName', 'Simulated')
        line(out.t./3600, ocv(out.za), 'Color', 'b', 'LineStyle', '--', 'DisplayName', 'Simulated')
        line(out.t./3600, ocv(out.zb), 'Color', 'r', 'LineStyle', '--', 'DisplayName', 'Simulated')
        xlabel('Time (hrs)')
        ylabel('Voltage (V)')
        legend show
    
        linkaxes([ax3, ax4, ax5, ax6], 'x')
        xlim(ax3, [0 t_max_hours])

    end

    fprintf('\n\n')

end
