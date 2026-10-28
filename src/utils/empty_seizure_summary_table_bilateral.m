function T = empty_seizure_summary_table_bilateral(tz)
% 0-row template for a RECONCILED seizures_summary.csv: empty_seizure_summary_table.m
% plus bilateral_columns('seizures_summary'). Used only when
% cfg.bilateral.rescue_mode ~= 'off' (see bilateral_reconcile.m).
    T = add_bilateral_columns(empty_seizure_summary_table(tz), 'seizures_summary');
end
