function [case_spec, note] = resolve_case_excel(attenuation_value, cfg)
% RESOLVE_CASE_EXCEL  Case of one channel from the recording Excel's
% Attenuation column (cfg.cases.from_excel = true), instead of the cases
% CSV / force / default that resolve_case.m uses.
%
%   [case_spec, note] = resolve_case_excel(attenuation_value, cfg)
%
% Attenuated (see excel_attenuation.m) -> profile
% cfg.cases.excel.case_attenuated, else cfg.cases.excel.case_clean.
% case_spec has the same fields as resolve_case.m's, with
% case_source = 'excel'. note is '' or a warning for qc_report.csv when
% the value was not in either list (treated as attenuated).

    [is_att, known] = excel_attenuation(attenuation_value, cfg);
    if is_att
        name = cfg.cases.excel.case_attenuated;
    else
        name = cfg.cases.excel.case_clean;
    end
    if ~isfield(cfg.cases.profiles, name)
        error('resolve_case_excel:BadCaseName', ...
            'cfg.cases.excel case "%s" is not defined in cfg.cases.profiles. Valid cases: %s', ...
            name, strjoin(fieldnames(cfg.cases.profiles), ', '));
    end
    profile = cfg.cases.profiles.(name);

    case_spec = struct('gain_mode', profile.gain_mode, 'gain', NaN, 'notch_mode', profile.notch_mode, ...
        'seizure_mode', profile.seizure_mode, 'case_applied', name, 'case_source', 'excel', ...
        'suggested_case', '');

    if known
        note = '';
    else
        note = sprintf('Excel Attenuation value "%s" not recognized: treated as attenuated', char(attenuation_value));
    end
end
