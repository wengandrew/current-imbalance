function elem(x, y, varargin)
    % Construct an elementary line object based on the 'line' function
    % This is just a wrapper for 'line' with a bunch of shortcuts to make
    % property-value pairs faster to instantiate

    p = inputParser;

    addOptional(p, 'ls', '-');
    addOptional(p, 'lw', 2);
    addOptional(p, 'm', 'none');
    addOptional(p, 'c', 'k');
    addOptional(p, 'ms', 5);
    addOptional(p, 'dn', '');
    addOptional(p, 'hv', 'on');
    addOptional(p, 'ax', gca);

    parse(p, varargin{:});

    r = p.Results;

    line(x, y, 'LineWidth', r.lw, ...
               'LineStyle', r.ls, ...
               'Color', r.c, ...
               'Marker', r.m, ...
               'MarkerFaceColor', r.c, ...
               'MarkerSize', r.ms, ...
               'DisplayName', r.dn, ...
               'Parent', r.ax, ...
               'HandleVisibility', r.hv)


end