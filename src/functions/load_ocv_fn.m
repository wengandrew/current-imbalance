function ocv = load_ocv_fn(type)
    % Build an OCV function from a data file
    % 
    % Parameters
    % ---------
    %   type: 'lfp' or 'nmc'
    %
    % Outputs
    % ---------
    %   ocv: an anonymous function OCV = f(z)

    switch type
        case 'lfp'
            filepath = 'data/ocv_Prada2013.csv';
        case 'nmc'
            filepath = 'data/ocv_Chen2020.csv';
        otherwise
            error('Type "%s" not supported.', type)
    end

    tbl = readtable(filepath);

    tbl.t = tbl.t - tbl.t(1);
    
    data.soc = tbl.t ./ max(tbl.t);
    data.ocv = tbl.V;

    ocv = @(z) interp1(data.soc, data.ocv, z, 'linear', 'extrap');

end