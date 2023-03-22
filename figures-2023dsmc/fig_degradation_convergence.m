function fig_degradation_convergence()
    % Run lifetime simulation for two parallel-connected batteries
    %
    % 2022/07/29
    % - Use the imbalance dynamics equations (inter-cycle dynamics)
    %   to simulate intra-cycle dynamics (degradation dynamics)
    % - At the end of each cycle, compute aggregate metrics and use it to
    %   update the capacities and resistances.
    % - Repeat.


    % Initialize the system
    p = initialize_cell_props_dsmc();

    [cyc_num_vec1, p1] = simulate(p, 0.5, 4);
    [cyc_num_vec2, p2] = simulate(p, 1.0, 0.002);
    [cyc_num_vec3, p3] = simulate(p, 1.5, 0.000002);


    set_default_plot_settings()

    fh = figure('Position', [500 100 1200 700]);
    th = tiledlayout(2, 3, 'Padding', 'none', 'TileSpacing', 'tight');
    th.TileIndexing = 'columnmajor';

    % Capacity loss
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);
    line(cyc_num_vec1, (p1.Qa - p1.dQa_vec)./3600, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec1, (p1.Qa - p1.dQa_control_vec)./3600, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', '--');
    line(cyc_num_vec1, (p1.Qb - p1.dQb_vec)./3600, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec1, (p1.Qb - p1.dQb_control_vec)./3600, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', '--');
    ylabel('$Q$ (Ah)', 'Interpreter', 'Latex')
    ylim([0 5])
    lh = legend('show'); set(lh, 'Location', 'NorthEast')

    ax2 = nexttile(th, 2); box on;
    line(cyc_num_vec1, (p1.Ra + p1.dRa_vec).*1000, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec1, (p1.Ra + p1.dRa_control_vec).*1000, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', '--');
    line(cyc_num_vec1, (p1.Rb + p1.dRb_vec).*1000, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec1, (p1.Rb + p1.dRb_control_vec).*1000, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', '--');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('$R$ (m$\Omega$)', 'Interpreter', 'Latex')
    lh = legend('show'); set(lh, 'Location', 'SouthEast');

    % Capacity loss
    ax3 = nexttile(th, 3); box on; set(ax3,'XTickLabel',[]);
    line(cyc_num_vec2, (p2.Qa - p2.dQa_vec)./3600, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec2, (p2.Qa - p2.dQa_control_vec)./3600, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', '--');
    line(cyc_num_vec2, (p2.Qb - p2.dQb_vec)./3600, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec2, (p2.Qb - p2.dQb_control_vec)./3600, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', '--');
    ylabel('$Q$ (Ah)', 'Interpreter', 'Latex')
    ylim([0 5])
    lh = legend('show'); set(lh, 'Location', 'NorthEast')

    ax4 = nexttile(th, 4); box on;
    line(cyc_num_vec2, (p2.Ra + p2.dRa_vec).*1000, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec2, (p2.Ra + p2.dRa_control_vec).*1000, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', '--');
    line(cyc_num_vec2, (p2.Rb + p2.dRb_vec).*1000, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec2, (p2.Rb + p2.dRb_control_vec).*1000, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', '--');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('$R$ (m$\Omega$)', 'Interpreter', 'Latex')
    lh = legend('show'); set(lh, 'Location', 'SouthEast');

    % Capacity loss
    ax5 = nexttile(th, 5); box on; set(ax5,'XTickLabel',[]);
    line(cyc_num_vec3, (p3.Qa - p3.dQa_vec)./3600, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec3, (p3.Qa - p3.dQa_control_vec)./3600, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', '--');
    line(cyc_num_vec3, (p3.Qb - p3.dQb_vec)./3600, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec3, (p3.Qb - p3.dQb_control_vec)./3600, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', '--');
    ylabel('$Q$ (Ah)', 'Interpreter', 'Latex')
    ylim([0 5])
    lh = legend('show'); set(lh, 'Location', 'NorthEast')

    ax6 = nexttile(th, 6); box on;
    line(cyc_num_vec3, (p3.Ra + p3.dRa_vec).*1000, 'Color', 'r', 'Marker', 'o', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
    line(cyc_num_vec3, (p3.Ra + p3.dRa_control_vec).*1000, 'Color', 'r', 'Marker', 'none', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A (control)', 'LineStyle', '--');
    line(cyc_num_vec3, (p3.Rb + p3.dRb_vec).*1000, 'Color', 'b', 'Marker', 'o', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
    line(cyc_num_vec3, (p3.Rb + p3.dRb_control_vec).*1000, 'Color', 'b', 'Marker', 'none', ...
        'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B (control)', 'LineStyle', '--');
    xlabel('Cycle Number', 'Interpreter', 'Latex')
    ylabel('$R$ (m$\Omega$)', 'Interpreter', 'Latex')
    lh = legend('show'); set(lh, 'Location', 'SouthEast');

%     subplot(223)
%     line(cyc_num_vec, Issa_vec, 'Color', 'r', 'Marker', 'o', ...
%         'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
%     line(cyc_num_vec, Issb_vec, 'Color', 'b', 'Marker', 'o', ...
%         'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
%     xlabel('Cycle Number', 'Interpreter', 'Latex')
%     ylabel('$I_{ss}$ (A)', 'Interpreter', 'Latex')
%     lh = legend('show');
%
%     subplot(224)
%     line(cyc_num_vec, zssa_vec, 'Color', 'r', 'Marker', 'o', ...
%         'MarkerFaceColor', 'r', 'MarkerSize', 2, 'DisplayName', 'Cell A');
%     line(cyc_num_vec, zssb_vec, 'Color', 'b', 'Marker', 'o', ...
%         'MarkerFaceColor', 'b', 'MarkerSize', 2, 'DisplayName', 'Cell B');
%     xlabel('Cycle Number', 'Interpreter', 'Latex')
%     ylabel('$z_{ss}$', 'Interpreter', 'Latex')
%     lh = legend('show');

end

function [cyc_num_vec, p] = simulate(p, n, gamma)

    cyc_num_vec = 0;

    % Define end condition
    dQa_vec = 0;
    dQb_vec = 0;
    dRa_vec = 0;
    dRb_vec = 0;
    dQa_control_vec = 0;
    dQb_control_vec = 0;
    dRa_control_vec = 0;
    dRb_control_vec = 0;
    Issa_vec = nan;
    Issb_vec = nan;
    zssa_vec = nan;
    zssb_vec = nan;

    % Initialize loop variables
    preva.dQ = 0;
    prevb.dQ = 0;
    preva.dR = 0;
    prevb.dR = 0;

    % Loop variables for the control cells (no current sharing)
    preva_control.dQ = 0;
    prevb_control.dQ = 0;
    preva_control.dR = 0;
    prevb_control.dR = 0;

    MAX_CYCLE_NUMBER = 2000;
    cyc_num = 1;

    za0 = 0;
    zb0 = 0;
    alpha = 1.2;
    U0 = 3.0;
    Vmin = 3.0;
    Vmax = 4.2;


%     f_ocv = @(z) U0 + alpha * z;
    f_ocv = load_ocv_fn('nmc');

    hours_to_charge = 3;
    current_target = -p.Qb / (hours_to_charge * 3600);
    t = linspace(0, 10*3600, 1.0e6)';
    I_chg = +current_target*ones(size(t)); % applied current (A)
    I_dch = -current_target*ones(size(t));
    I_current_cutoff = current_target/20;

    % Loop over cycles
    while cyc_num < MAX_CYCLE_NUMBER

        fprintf('Simulating cycle %g ...\n', cyc_num)

        % Solve for the intra-cycle dynamics
        %
        % curr: struct holding current (Ra, Rb, Qa, Rb)
        % prev: struct holding previous (Ra, Rb, Qa, Rb)
        Qa = p.Qa - preva.dQ;
        Qb = p.Qb - prevb.dQ;
        Ra = p.Ra + preva.dR;
        Rb = p.Rb + prevb.dR;

        % Update capacity loss and resistance growth increment for this cycle
        delta_time = hours_to_charge * 3600;

        [Issa, Issb, zssa, zssb] = update_cycle_metrics(Qa, Qb, Ra, Rb, alpha, ...
            I_chg, I_dch, I_current_cutoff, za0, zb0, delta_time, ...
            f_ocv, Vmin, Vmax);

        curra = update_states_iss(preva, delta_time, Issa, n, gamma);
        currb = update_states_iss(prevb, delta_time, Issb, n, gamma);

        curra_control = update_states_iss(preva_control, delta_time, I_chg(1)/2, n, gamma);
        currb_control = update_states_iss(prevb_control, delta_time, I_chg(1)/2, n, gamma);

        dQa_vec = [dQa_vec ; curra.dQ];
        dQb_vec = [dQb_vec ; currb.dQ];
        dRa_vec = [dRa_vec ; curra.dR];
        dRb_vec = [dRb_vec ; currb.dR];

        dQa_control_vec = [dQa_control_vec ; curra_control.dQ];
        dQb_control_vec = [dQb_control_vec ; currb_control.dQ];
        dRa_control_vec = [dRa_control_vec ; curra_control.dR];
        dRb_control_vec = [dRb_control_vec ; currb_control.dR];

        Issa_vec = [Issa_vec ; Issa];
        Issb_vec = [Issb_vec ; Issb];
        zssa_vec = [zssa_vec ; zssa];
        zssb_vec = [zssb_vec ; zssb];

        preva = curra;
        prevb = currb;
        preva_control = curra_control;
        prevb_control = currb_control;

        cyc_num_vec = [cyc_num_vec; cyc_num];
        cyc_num = cyc_num + 1;

        if Qa < 0.0 * 3600
            break
        end
        if Qb < 0.0 * 3600
            break
        end

    end

    p.dQa_vec = dQa_vec;
    p.dQb_vec = dQb_vec;
    p.dRa_vec = dRa_vec;
    p.dRb_vec = dRb_vec;

    p.issa_vec = Issa_vec;
    p.issb_vec = Issb_vec;
    p.zssa_vec = zssa_vec;
    p.zssb_vec = zssb_vec;

    p.dQa_control_vec = dQa_control_vec;
    p.dQb_control_vec = dQb_control_vec;
    p.dRa_control_vec = dRa_control_vec;
    p.dRb_control_vec = dRb_control_vec;

end

function [Issa, Issb, zssa, zssb] = update_cycle_metrics(Qa, Qb, Ra, Rb, ...
            alpha, I_chg, I_dch, I_current_cutoff, za0, zb0, dt, f_ocv, Vmin, Vmax)
    % Compute steady-state values

    type = 'analytic';

    switch type
        case 'analytic'
            Q = Qa + Qb;
            R = Ra + Rb;

            kappa = 1/alpha * (Ra*Qa - Rb*Qb) / Q;
            tau = R/alpha * (Qa*Qb / (Qa+Qb));
            Issa = 1/R * ((Rb + alpha*kappa))*I_chg(1);
            Issb = 1/R * ((Ra - alpha*kappa))*I_chg(1);

            zssa = 1/Q * (Qa*za0 + Qb*zb0) + 1/Q * (+Qb*kappa*(1-exp(-dt/tau)) - dt) * I_chg(1);
            zssb = 1/Q * (Qa*za0 + Qb*zb0) + 1/Q * (-Qa*kappa*(1-exp(-dt/tau)) - dt) * I_chg(1);

            fprintf('%.3fA, %.3fA, %.3f, %.3f\n', Issa, Issb, zssa, zssb)

        case 'simulation'
            % Run a discrete-time simulation with a non-linear OCV function
            out = run_discrete_time_simulation_complete(I_chg, I_dch, ...
               I_current_cutoff, Qa, Qb, Ra, Rb, za0, zb0, f_ocv, ...
               Vmin, Vmax);

            % Here, 'steady-state' is more like 'final value', since a true
            % steady-state doesn't exist when the OCV function is changing
            % all over the place due to the non-linearities in the
            % function. The 'final value' is the final value at the end of
            % CC phase only, since this point represents the state at which
            % the cell is most likely to lithium plate.

            idx_cc = find(out.state == 0);
            Issa = out.Ia(idx_cc); Issa = Issa(end);
            Issb = out.Ib(idx_cc); Issb = Issb(end);
            zssa = out.za(idx_cc); zssa = zssa(end);
            zssb = out.zb(idx_cc); zssb = zssb(end);

            fprintf('Sim | %.3fA, %.3fA, %.3f, %.3f\n', Issa, Issb, zssa, zssb)

    end

%             Q = Qa + Qb;
%             R = Ra + Rb;
%
%             kappa = 1/alpha * (Ra*Qa - Rb*Qb) / Q;
%
%             Issa = 1/R * ((Rb + alpha*kappa))*I_chg(1);
%             Issb = 1/R * ((Ra - alpha*kappa))*I_chg(1);
%
%             zssa = 1/Q * (Qa*za0 + Qb*zb0) + 1/Q * (+Qb*kappa - dt) * I_chg(1);
%             zssb = 1/Q * (Qa*za0 + Qb*zb0) + 1/Q * (-Qa*kappa - dt) * I_chg(1);
%
%             fprintf('Ana | %.3fA, %.3fA, %.3f, %.3f\n', Issa, Issb, zssa, zssb)

end


function curr = update_states_iss(prev, delta_time, Iss, n, gamma)
    % Q-R state update function based on steady-state current imbalance
    % values
    %
    % Parameters
    % ---------
    % prev (struct)
    %    previous state (includes dQ and dR)
    % delta_time (float)
    %    time for the cycle
    % Iss (float)
    %    steady-state current value in Amperes
    %
    % Outputs
    % ---------
    % curr
    %    updated state

    % Update rate constants for capacities and resistances
    kq = gamma*abs(Iss);

    kr = 0.000002*kq;

    % Discrete state update equation for SEI growth
    curr.dQ = ( kq^(1/n) * delta_time + prev.dQ^(1/n) )^n;
    curr.dR = ( kr^(1/n) * delta_time + prev.dR^(1/n) )^n;

end
