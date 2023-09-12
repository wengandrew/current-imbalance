function aw20221221_imbalance_bound_simulations(analysis_type, q, r, current_target)
% 
% 
% Args:
%  analysis_type : 'lfp', 'nmc', 'quad'ratic function
%  q : Qa/Qb
%  r : Ra/Rb
%  current_target : target current in amperes

    to_plot = true;

    set_default_plot_settings()

    switch analysis_type
        case 'lfp'
            alpha = 0.6;
            Vmax = 3.6;
            Vmin = 3.0;
            ocv = load_ocv_fn('lfp');
        case 'nmc'
            alpha = 1.2;
            Vmax = 4.2;
            Vmin = 3.0;
            ocv = load_ocv_fn('nmc');
        case 'quad'
            alpha = 1.2;
            Vmax = 4.2;
            Vmin = 3.0;
            ocv = @(x) 0.6.*x.^2 + 0.6.*x + 3;
    end

    % Initialize model parameters
    Ra = 0.15; % ohms
    Qa = 3 * 3600; % As
    za0 = 0.00;
    zb0 = 0.00;

    dz0 = za0 - zb0;

    Qb = Qa/q;
    Rb = Ra/r;

    p.Ra = Ra;
    p.Qa = Qa;
    p.Rb = Rb;
    p.Qb = Qb;
        
    kappa = 1 / alpha * (Ra*Qa - Rb*Qb) / (Qa+Qb);

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
        alpha, Ra, Rb, Qa, Qb, za0, zb0, Vmax, Vmin);

    %% Test the simulation based solution
    out = run_discrete_time_simulation_cccv(t, I, Qa, Qb, Ra, Rb, ...
            za0, zb0, ocv, Vmax);

    % SOC imbalance bounds
    [condition, zbound_l2, zbound_linf, ibound_l2, ibound_linf] = ...
        solve_imbalance_bounds(ocv, p, out.t, I, dz0);

    z_tilde = out.zb - out.za;
    signorm_z_tilde = sqrt(trapz(out.t, z_tilde.^2));

%     fprintf('L2   Condition:  %.3f <= %.3f ? \t%g \n', abs(current_target), condition, abs(current_target) < condition);
%     fprintf('L2   Bound:      %.3f <= %.3f ? \t%g \n', signorm_z_tilde, zbound_l2, signorm_z_tilde < zbound_l2)
    
%     fprintf('Linf Condition:  %.3f <= %.3f ? \t%g \n', max(abs(I)), condition, max(abs(I)) < condition);
%     fprintf('Linf Bound:      %.3f <= %.3f ? \t%g \n', max(dz0), zbound_linf, max(dz0) < zbound_linf);

    % See the results

    if to_plot

        figure('OuterPosition', [0, 0, 900, 1200])

        ax3 = subplot(511);
        line(tfinal./3600, za, 'Color', 'b', 'DisplayName', 'Analytic')
        line(tfinal./3600, zb, 'Color', 'r', 'DisplayName', 'Analytic')
        line(out.t./3600, out.za, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
        line(out.t./3600, out.zb, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
        ylabel('$z$', 'Interpreter', 'Latex')
        legend showline


        ax4 = subplot(512);
        line(tfinal./3600, abs(za - zb), 'Color', 'r', 'DisplayName', 'Analytic')
        yline(abs(kappa * current_target), 'DisplayName', 'Bound ($\kappa I$)', 'LineStyle', '--', 'Color', 'r')
        line(out.t./3600, abs(out.za - out.zb), 'Color', 'b', 'DisplayName', 'Simulated')
        line(out.t./3600, zbound_linf, 'DisplayName', 'Bound ($-\frac{B}{Ak_1}\max(|I|) + |\tilde{z}_0|$)', 'LineStyle', '--', 'Color', 'b')
        ylabel('$|\Delta z|$', 'Interpreter', 'Latex')
        legend show
    
        ax5 = subplot(513);
        line(tfinal./3600, Ia, 'Color', 'b', 'DisplayName', 'Analytic')
        line(tfinal./3600, Ib, 'Color', 'r', 'DisplayName', 'Analytic')
        line(out.t./3600, out.Ia, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Simulated')
        line(out.t./3600, out.Ib, 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Simulated')
        ylabel('Current (A)')
        legend show

        ax5 = subplot(514);
        line(tfinal./3600, abs(Ia-Ib), 'Color', 'r', 'DisplayName', 'Analytic')
        yline(abs((Qa-Qb)/(Qa+Qb)*current_target), 'Color', 'r', 'LineStyle', ':', 'DisplayName', 'Bound ($\frac{Q_a-Q_b}{Q_a+Q_b}I$)')
        line(out.t./3600, out.Ia - out.Ib, 'Color', 'b', 'DisplayName', 'Simulated')
        yline(ibound_linf, 'Color', 'b', 'LineStyle', ':', 'DisplayName', 'Bound')
        ylabel('$|\Delta I|$ (A)', 'Interpreter', 'Latex')
        legend show
    
        ax6 = subplot(515);
        line(tfinal/3600, Vmin + alpha*za - Ia*Ra, 'Color', 'k', 'DisplayName', 'Analytic')
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
