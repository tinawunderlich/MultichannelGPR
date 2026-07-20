function data = readGroup(radarFile)
    % reads complete group of radargrams_*.mat/global_coords_*.mat and
    % x_*.mat (give path to radargrams_*.mat)

    % get index from filename
    [folder, name, ~] = fileparts(radarFile);
    idx = regexp(name, '\d+$', 'match', 'once');  % z.B. "1", "2", ...

    % load radargrams and reconstruct cell array if necessary
    r = load(radarFile);
    varnames = fieldnames(r);
    if isscalar(varnames) && strcmp(varnames{1},'radargrams')
        radargrams=r.radargrams;
    else
       varnames = sort(fieldnames(r.rg));
        radargrams = cell(length(varnames), 1);
        for i = 1:length(varnames)
            radargrams{i} = r.rg.(varnames{i});
        end
    end
    c = load(fullfile(folder, ['global_coords_' idx '.mat']));
    x = load(fullfile(folder, ['x_' idx '.mat']));

    % give back as struct
    data.radargrams    = radargrams;
    data.global_coords = c.global_coords;
    data.x             = x.x;
    data.index         = str2double(idx);
end