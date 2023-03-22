function ocv = load_ocv_fn(type)
    % Build an OCV function from a data file
    % 
    % Parameters
    % ---------
    %   type: 'lfp' or 'nmc' or 'nmc-umbl2022feb'
    %
    % Outputs
    % ---------
    %   ocv: an anonymous function OCV = f(z)

    switch type
        case 'lfp'
            filepath = 'data/ocv_Prada2013.csv';
        case 'nmc'
            filepath = 'data/ocv_Chen2020.csv';
        case 'nmc-umbl2022feb'
            filepath = 'data/20230303_c20_charge.csv';
        otherwise
            error('Type "%s" not supported.', type)
    end

    tbl = readtable(filepath);

    tbl.t = tbl.t - tbl.t(1);
    
    data.soc = tbl.t ./ max(tbl.t);
    data.ocv = tbl.V;

    ocv = griddedInterpolant(data.soc, data.ocv, 'linear', 'nearest');

end
