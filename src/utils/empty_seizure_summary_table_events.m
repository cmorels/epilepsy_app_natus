function T = empty_seizure_summary_table_events(tz, with_events)
% 0-row template for seizures_summary.csv with the review band on:
% empty_seizure_summary_table.m, plus event_columns('seizures_summary') when
% the run had robust channels (with_events). The band-off template is left
% unchanged so band-off runs stay byte-identical.
    T = empty_seizure_summary_table(tz);
    if with_events
        T = add_event_columns(T, 'seizures_summary');
    end
end
