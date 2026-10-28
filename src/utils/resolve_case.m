function [case_spec, matched_row_idx] = resolve_case(candidate_names, region, cases_table, cfg)
% RESOLVE_CASE  Resolve the case (gain/notch mode) for one channel, by
% priority: (1) explicit gain/notch values on its matching CSV row,
% (2) that row's case name, (3) cfg.cases.force, (4) cfg.cases.default.
% Explicit gain/notch values override their OWN field only -- e.g. a row
% with case='attenuated' (gain=auto, notch=off) and an explicit
% notch='on' resolves to gain_mode='auto' (from the case), notch_mode='on'
% (explicit override), case_applied='attenuated', case_source='explicit'.
%
%   [case_spec, case_applied, case_source, matched_row_idx] = ...
%       resolve_case(candidate_names, region, cases_table, cfg)
%
%   candidate_names : cellstr of names this channel may be referred to by
%                      in the CSV's source_file column (typically the EDF
%                      filename and this channel's own txt filename --
%                      "source_file admite el nombre del EDF o el del
%                      txt")
%   region          : this channel's region (matched against the row's
%                      region column; a row with an empty region column
%                      matches every region in that file)
%   cases_table     : table from load_cases.m, or [] / 0-row if no CSV
%                      was configured (cfg.cases.file == '')
%   cfg             : struct from pipeline_config.m
%
% OUTPUTS: case_spec (struct: gain_mode, gain, notch_mode, case_applied --
% the resolved case NAME -- and case_source, one of 'explicit'|'csv'|
% 'force'|'default'), matched_row_idx (row index into cases_table that
% matched, or NaN if none -- used by the caller to track which CSV rows
% were consumed, for the "row never matched a file" warning).

    row = [];
    matched_row_idx = NaN;
    if ~isempty(cases_table) && height(cases_table) > 0
        [row, matched_row_idx] = find_matching_row(cases_table, candidate_names, region);
    end

    if ~isempty(row)
        row_case = strtrim(row.('case'){1});
        row_gain = row.gain(1);
        row_notch = strtrim(row.notch{1});
        if ismember('suggested_case', row.Properties.VariableNames)
            suggested_case = row.suggested_case{1};
        else
            suggested_case = '';
        end
    else
        row_case = '';
        row_gain = NaN;
        row_notch = '';
        suggested_case = '';
    end

    if ~isempty(row_case)
        base_name = row_case;
        base_source = 'csv';
    elseif ~isempty(cfg.cases.force)
        base_name = cfg.cases.force;
        base_source = 'force';
    else
        base_name = cfg.cases.default;
        base_source = 'default';
    end

    if ~isfield(cfg.cases.profiles, base_name)
        error('resolve_case:BadCaseName', ...
            'Case "%s" is not defined in cfg.cases.profiles. Valid cases: %s', ...
            base_name, strjoin(fieldnames(cfg.cases.profiles), ', '));
    end
    profile = cfg.cases.profiles.(base_name);

    gain_mode = profile.gain_mode;
    gain_value = NaN;
    notch_mode = profile.notch_mode;
    seizure_mode = profile.seizure_mode;
    explicit_used = false;

    if ~isnan(row_gain)
        gain_mode = 'explicit';
        gain_value = row_gain;
        explicit_used = true;
    end
    if ~isempty(row_notch)
        nm = lower(row_notch);
        if ~ismember(nm, {'on', 'off'})
            error('resolve_case:BadNotchValue', ...
                'notch column must be "on" or "off" (got "%s") for a row matching %s.', ...
                row_notch, strjoin(candidate_names, ' / '));
        end
        notch_mode = nm;
        explicit_used = true;
    end

    if explicit_used
        case_source = 'explicit';
    else
        case_source = base_source;
    end
    case_applied = base_name;

    case_spec = struct('gain_mode', gain_mode, 'gain', gain_value, 'notch_mode', notch_mode, ...
        'seizure_mode', seizure_mode, 'case_applied', case_applied, 'case_source', case_source, ...
        'suggested_case', suggested_case);
end

function [row, idx] = find_matching_row(cases_table, candidate_names, region)
    file_mask = false(height(cases_table), 1);
    for i = 1:numel(candidate_names)
        file_mask = file_mask | strcmpi(cases_table.source_file, candidate_names{i});
    end

    region_specific = file_mask & strcmpi(cases_table.region, region) & ~cellfun(@isempty, cases_table.region);
    region_wildcard = file_mask & cellfun(@isempty, cases_table.region);

    if any(region_specific)
        candidates = find(region_specific);
    elseif any(region_wildcard)
        candidates = find(region_wildcard);
    else
        row = [];
        idx = NaN;
        return;
    end

    idx = candidates(end);  % last matching row wins if more than one
    row = cases_table(idx, :);
end
