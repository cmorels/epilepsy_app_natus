function T = empty_seizure_events_table_events(tz, with_events)
% 0-row template for seizures_events.csv with the review band on:
% empty_seizure_events_table.m + ll_status, plus event_columns('seizures_events')
% when the run had robust channels (with_events). The band-off template
% (empty_seizure_events_table.m) is deliberately left unchanged so band-off
% runs stay byte-identical.
    T = add_event_columns(empty_seizure_events_table(tz), 'll_status');
    if with_events
        T = add_event_columns(T, 'seizures_events');
    end
end
