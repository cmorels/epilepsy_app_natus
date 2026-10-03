function T = empty_seizure_events_table_bilateral(tz)
% 0-row template for a RECONCILED seizures_events.csv: empty_seizure_events_table.m
% plus bilateral_columns('seizures_events'). Used only when
% cfg.bilateral.rescue_mode ~= 'off' (see bilateral_reconcile.m).
    T = add_bilateral_columns(empty_seizure_events_table(tz), 'seizures_events');
end
