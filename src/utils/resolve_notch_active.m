function notch_active = resolve_notch_active(notch_mode, quality_class)
% RESOLVE_NOTCH_ACTIVE  Whether the notch actually runs for this channel,
% given its resolved notch_mode and (for 'auto') its quality_class. Shared
% by clean_lfp.m (which applies it) and run_pipeline_edf.m's selective-
% reprocessing check (which needs to predict it WITHOUT re-running
% clean_lfp, to decide whether re-running is even necessary) -- kept in
% one place so the two can't drift apart.
    switch notch_mode
        case 'off'
            notch_active = false;
        case 'on'
            notch_active = true;
        case 'auto'
            notch_active = ismember(quality_class, {'line_contaminated', 'low_amplitude+line'});
        otherwise
            error('resolve_notch_active:BadMode', ...
                'notch_mode must be ''off'', ''on'', or ''auto'' (got ''%s'').', notch_mode);
    end
end
