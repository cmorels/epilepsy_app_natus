function [is_attenuated, known] = excel_attenuation(value, cfg)
% EXCEL_ATTENUATION  Translate one Attenuation cell of the recording Excel
% (EEG_recording_log.xlsx) into attenuated yes/no.
%
%   [is_attenuated, known] = excel_attenuation(value, cfg)
%
% value is compared lower-cased, trimmed and with runs of spaces collapsed
% against cfg.cases.excel.attenuation_no (default '', 'no', 'none') and
% cfg.cases.excel.attenuation_yes (default 'severe', 'mild',
% 'unsure, perhaps mild'). A value in neither list is treated as
% attenuated (known = false) so the caller can warn about it -- never an
% error.

    v = lower(strtrim(regexprep(char(value), '\s+', ' ')));
    no_list = lower(strtrim(cfg.cases.excel.attenuation_no));
    yes_list = lower(strtrim(cfg.cases.excel.attenuation_yes));
    if ismember(v, no_list)
        is_attenuated = false;
        known = true;
    elseif ismember(v, yes_list)
        is_attenuated = true;
        known = true;
    else
        is_attenuated = true;
        known = false;
    end
end
