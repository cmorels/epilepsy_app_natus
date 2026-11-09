function o = rp128_paths()
% RP128_PATHS  Paths and constants of the reprocessing of 128-s, EDF
% 20260729~ TRAP_d5154696 (C2): 128-s moved from port A1 to A2 during the
% recording; A1 is valid until 17:18:39, A2 from 17:23:55 (gap = DX-21).
% Used by every rp128_* function.

    o.tz = 'Europe/Paris';
    o.batch = '\\iws\epic\users\camila.morel\2-Analysis\ECoG\batch 2026';
    o.repo = fullfile(o.batch, 'epilepsy_app_natus');
    o.edf_dir = '\\iws\epic\users\camila.morel\1-DATA\Video_ECoG\2026\natus\EDF export\decimated';
    o.subject = '128-s';
    o.stem = '20260729~ TRAP_d5154696-610f-41eb-91af-d50e5658bcd3';
    o.src_edf = fullfile(o.edf_dir, o.subject, [o.stem '.EDF']);            % new export: A1 and A2, whole recording
    o.next_stem = '20260730~ TRAP_02be1c26-73a5-4363-93fe-ac7f6e5416e8';
    o.next_edf = fullfile(o.edf_dir, o.subject, [o.next_stem '.EDF']);
    o.master_log = fullfile(o.edf_dir, 'EEG_recording_log.xlsx');
    o.c2_root = fullfile(o.batch, 'pipeline_output_campaign_20261014');
    o.c2_subject = fullfile(o.c2_root, o.subject);
    o.c2_txt = struct('HPCr', fullfile(o.c2_subject, '01_txt', '128-s_20260729_100751_HPCr.txt'), ...
        'HPCl', fullfile(o.c2_subject, '01_txt', '128-s_20260729_100751_HPCl.txt'));   % original EDF (A1) as the C2 imported it
    o.c2_ref_csv = fullfile(o.c2_root, 'reference_per_region.csv');
    o.c2_cfg_mat = fullfile(o.c2_subject, 'config_used.mat');

    o.out = fullfile(o.batch, o.subject, 'reprocess_20260729_A1A2');
    o.log_override = fullfile(o.out, 'EEG_recording_log_override.xlsx');
    % decided 2026-10-30 (option 1): two EDFs, A1 only and A2 only, so each
    % gets its own detection thresholds (A2 has ~2x A1's broadband amplitude)
    o.seg = struct('name', {'A1', 'A2'}, 'stem', {[o.stem '_128s_A1'], [o.stem '_128s_A2']}, ...
        'port', {'A1', 'A2'}, 'labels', {{'EEG A1C2', 'EEG A1C4'}, {'EEG A2C2', 'EEG A2C4'}});
    for i = 1:2, o.seg(i).edf = fullfile(o.out, [o.seg(i).stem '.EDF']); end
    o.spr = 512;        % 1 s records: the cut falls at most 1 s before each segment's end
    o.ann_spr = 128;    % 256 bytes of annotations per record

    o.cut1 = datetime(2026, 7, 29, 17, 18, 39, 'TimeZone', o.tz);   % DX-21 start: A1 strictly before
    o.cut2 = datetime(2026, 7, 29, 17, 23, 55, 'TimeZone', o.tz);   % DX-21 end: A2 from here on
    o.ports = {'A1', 'A2'};
    o.regions = {'HPCr', 'HPCl'};
    o.compat_window_s = 3600;
    o.compat_ratio = [0.8 1.25];
end
