function fig_degradation_sei()
    % Run lifetime simulation with the discrete SEI growth model
    %
    % 2022/08/03
    % - just demonstrate how we can use the discrete update version of the
    %   SEI growth model to handle changing use cases.

    % Initialize the system
    p.Qa = 3 * 3600; % Amp-seconds
    p.Ra = 150 / 1000; % Ohms

    p.U0 = 3.0;
    p.alpha = 1.2;

    p.f_Ua = @(z) U0 + alpha * z;
    p.f_Ub = @(z) U0 + alpha * z;

    % Define end condition
    dQa_vec = 0;
    dRa_vec = 0;
    dQa1_vec = 0;
    dRa1_vec = 0;
    dQa2_vec = 0;
    dRa2_vec = 0;
    kc_vec = 0;

    % Initialize loop variables
    preva.dQ = 0;
    preva.dR = 0;
    preva1.dQ = 0;
    preva1.dR = 0;
    preva2.dQ = 0;
    preva2.dR = 0;

    MAX_CYCLE_NUMBER = 1000;
    cyc_num_vec = 0;
    cyc_num = 1;

    hours_to_charge = 3;

    % SEI reaction rate constants
    ka = 1;
    kb = 2;

    % Loop over cycles
    while cyc_num < MAX_CYCLE_NUMBER

        fprintf('Simulating cycle %g ...\n', cyc_num)

        if cyc_num < 300
            kc = ka;
        elseif cyc_num > 600
            kc = ka;
        else
            kc = kb;
        end

        % Update capacity loss and resistance growth increment for this cycle
        delta_time = hours_to_charge * 3600;

        curra = update_states_simple(preva, delta_time, kc);
        curra1 = update_states_simple(preva1, delta_time, ka);
        curra2 = update_states_simple(preva2, delta_time, kb);

        kc_vec = [kc_vec ; kc];
        dQa_vec = [dQa_vec ; curra.dQ];
        dRa_vec = [dRa_vec ; curra.dR];
        dQa1_vec = [dQa1_vec ; curra1.dQ];
        dRa1_vec = [dRa1_vec ; curra1.dR];
        dQa2_vec = [dQa2_vec ; curra2.dQ];
        dRa2_vec = [dRa2_vec ; curra2.dR];
        preva = curra;
        preva1 = curra1;
        preva2 = curra2;

        cyc_num_vec = [cyc_num_vec; cyc_num];
        cyc_num = cyc_num + 1;

    end

    set_default_plot_settings()


    fh = figure('Position', [500 100 600 600]);
    th = tiledlayout(3, 1, 'Padding', 'none', 'TileSpacing', 'none');

    % Capacity loss
    ax1 = nexttile(th, 1); box on; set(ax1,'XTickLabel',[]);
    ylabel('$Q_i$ (Ah)', 'Interpreter', 'Latex')

    line(cyc_num_vec, (p.Qa - dQa_vec)./3600, 'Color', 'k',  ...
        'LineStyle', '-', 'LineWidth', 3, 'DisplayName', 'Mixed');
    line(cyc_num_vec, (p.Qa - dQa1_vec)./3600, 'Color', 'r', ...
        'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Rate 1');
    line(cyc_num_vec, (p.Qa - dQa2_vec)./3600, 'Color', 'b', ...
        'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Rate 2');
    ylim([0.8 3])
    lh = legend('show'); set(lh, 'Location', 'SouthWest')

    % Resistance growth
    ax2 = nexttile(th, 2); set(ax2, 'XTickLabel', []); box on;
    ylabel(ax2, '$R_i$ (m$\Omega$)', 'Interpreter', 'Latex')

    line(cyc_num_vec, (p.Ra + dRa_vec).*1000, 'Color', 'k', ...
        'LineStyle', '-', 'LineWidth', 3, 'DisplayName', 'Mixed');
    line(cyc_num_vec, (p.Ra + dRa1_vec).*1000, 'Color', 'r', ...
        'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Rate 1');
    line(cyc_num_vec, (p.Ra + dRa2_vec).*1000, 'Color', 'b', ...
        'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Rate 2');

    xlabel('Cycle Number', 'Interpreter', 'Latex')
%     lh = legend('show', 'Interpreter', 'Latex');

    % Rate constant
    ax3 = nexttile(th, 3); box on;
    ylabel(ax3, '$r_i$', 'Interpreter', 'Latex')
    xlabel(ax3, 'Cycle Number')

    line(cyc_num_vec, kc_vec, 'Color', 'k', ...
        'LineStyle', '-', 'LineWidth', 3, ...
        'DisplayName', 'Mixed');
    line(cyc_num_vec, ka*ones(size(cyc_num_vec)), 'Color', 'r', ...
        'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Rate 1');
    line(cyc_num_vec, kb*ones(size(cyc_num_vec)), 'Color', 'b', ...
        'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Rate 2');

    xlabel('$n$', 'Interpreter', 'Latex')
%     lh = legend('show', 'Interpreter', 'Latex');
    ylim([min([ka, kb, kc])*0.9 max([ka, kb, kc])*1.1])


end

function curr = update_states_simple(prev, delta_time, kc)
    % Q-R state update function
    %
    % Parameters
    % ---------
    % prev (struct)
    %    previous state (includes dQ and dR)
    % delta_time (float)
    %    time for the cycle
    % kc (float)
    %    reaction rate constant for capacity loss based on SEI growth
    %
    % Outputs
    % ---------
    % curr
    %    updated state

    n = 0.5;

    % Rate constant for resistance growth is proportional to rate constant for
    % capacity loss
    kr = kc./1e5;

    % Discrete state update equation for SEI growth
    curr.dQ = (kc^(1/n) * delta_time + prev.dQ^(1/n))^n;
    curr.dR = (kr^(1/n) * delta_time + prev.dR^(1/n))^n;

end
