function rp128_prepare_labels()
% RP128_PREPARE_LABELS  Phase 5: leave 05_summaries/seizures_events.csv of
% the reprocessing ready to label, like the campaigns' reviewed files:
% the pipeline's 30 columns + AE 'user_check' and AF 'comments', empty.
% The pipeline's own file is kept untouched as seizures_events_pipeline.csv
% (full precision; run_campaign_stats joins the labels to it row by row).
% Refuses to touch a seizures_events.csv that already has a user_check.

    o = rp128_paths();
    d = fullfile(o.out, '05_summaries');
    lab = fullfile(d, 'seizures_events.csv');
    pri = fullfile(d, 'seizures_events_pipeline.csv');
    L = readlines(lab, 'EmptyLineRule', 'skip');
    if contains(L(1), 'user_check')
        error('rp128_prepare_labels:Labeled', '%s already has a user_check column: not touched', lab);
    end
    if ~isfile(pri), copyfile(lab, pri); end
    P = readlines(pri, 'EmptyLineRule', 'skip');
    n = numel(split(P(1), ','));
    if n ~= 30, error('rp128_prepare_labels:Columns', '%s has %d columns, expected the pipeline''s 30', pri, n); end
    out = [P(1) + ",user_check,comments"; P(2:end) + ",,"];
    fid = fopen(lab, 'w');
    fprintf(fid, '%s\n', out);
    fclose(fid);
    fprintf('%s: %d detection row(s), columns AE user_check / AF comments empty\n', lab, numel(P) - 1);
end
