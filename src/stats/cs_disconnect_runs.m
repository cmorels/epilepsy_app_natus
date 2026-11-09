function runs = cs_disconnect_runs(evf, max_gap_s)
% CS_DISCONNECT_RUNS  Rule R5: FP events labeled 'disconnect' that follow
% each other (time order, no other event in between) AND are separated by
% at most max_gap_s form one zone from the first start to the last end.
% Without the gap bound, "consecutive" would join disconnect detections
% hours or days apart (e.g. 070-s 20260807: 2 events 49 h apart).
% Isolated disconnect events are returned as their own span.
%
%   evf  : events of ONE (campaign, mouse, EDF), any label
%   runs : Nx2 POSIX s

    evf = sortrows(evf, 'start_posix');
    isd = strcmp(evf.label, 'FP') & contains(evf.comments, 'disconnect');
    runs = zeros(0, 2);
    i = 1;
    while i <= height(evf)
        if ~isd(i)
            i = i + 1;
            continue;
        end
        j = i;
        while j < height(evf) && isd(j+1) && evf.start_posix(j+1) - max(evf.end_posix(i:j)) <= max_gap_s
            j = j + 1;
        end
        runs(end+1, :) = [evf.start_posix(i) max(evf.end_posix(i:j))]; %#ok<AGROW>
        i = j + 1;
    end
end
