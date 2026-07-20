function data = readRadargrams(radarFile)
    % reads radargrams_*.mat and optionally reconstructs cell structure

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

    % give back as struct
    data.radargrams    = radargrams;
    data.index         = str2double(idx);
end