function paths = make_spec_session(folder, subject, session_start, region, x, fs, varargin)
% MAKE_SPEC_SESSION  Test helper: writes a pipeline-shaped clean txt (and
% optionally its raw 01_txt) with write_lfp_txt.m.
%
%   paths = make_spec_session(folder, subject, session_start, region, x, fs, ...
%       'source_file', 'EDFA.EDF', 'raw_folder', '', 'gain', 1, 'n_samples_header', [])
%
% session_start: zoned datetime. The clean header carries the fields
% clean_lfp.m writes that the module reads (case/gain/notch/thresholds).
% With 'raw_folder', the raw txt holds x / gain (so raw * gain == clean).
% 'n_samples_header' lets a header claim a longer recording than x (only
% used by tests that never read the data).

    p = inputParser;
    p.addParameter('source_file', 'EDFA.EDF');
    p.addParameter('raw_folder', '');
    p.addParameter('gain', 1);
    p.addParameter('n_samples_header', []);
    p.parse(varargin{:});
    o = p.Results;

    if ~isfolder(folder)
        mkdir(folder);
    end
    ss = session_start;
    ss.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX';
    base = sprintf('%s_%s_%s', subject, char(string(session_start, 'yyyyMMdd')), char(string(session_start, 'HHmmss')));
    n = numel(x);
    nh = n;
    if ~isempty(o.n_samples_header)
        nh = o.n_samples_header;
    end
    hdr = { ...
        'subject_id', subject; 'fs', fs; 'time_unit', 'seconds'; 'region', region; ...
        'channel_label', ['EEG ' region]; 'source_file', o.source_file; 'source_format', 'EDF+C'; ...
        'session_start_datetime', char(ss); 'session_start_unix', posixtime(session_start); ...
        'timezone', session_start.TimeZone; 'n_samples', nh; 'n_valid_samples', sum(~isnan(x)); ...
        'n_gaps', double(any(isnan(x))); 'total_duration_s', nh / fs; 'units', 'microvolts'};
    paths = struct('clean', '', 'raw', '');
    if ~isempty(o.raw_folder)
        if ~isfolder(o.raw_folder)
            mkdir(o.raw_folder);
        end
        paths.raw = fullfile(o.raw_folder, [base '_' region '.txt']);
        write_lfp_txt(paths.raw, hdr, x / o.gain);
    end
    clean = [hdr; { ...
        'clean_source_file', [base '_' region '.txt']; 'threshold_type', 'FIXED (user config)'; ...
        'lower_threshold_uV', -2500; 'upper_threshold_uV', 2500; 'n_outliers', 0; 'outlier_pct', 0; ...
        'case_applied', 'excel_clean'; 'gain_applied', o.gain; 'gain_source', 'explicit'; ...
        'quality_class', 'ok'; 'notch_applied', false; 'notch_freqs', ''; 'notch_halfwidth_hz', 1; 'notch_order', 4}];
    paths.clean = fullfile(folder, [base '_' region '_clean.txt']);
    write_lfp_txt(paths.clean, clean, x);
end
