function T = apply_tz(T, tz)
% APPLY_TZ  Stamp tz onto every datetime column of T.
%
% All-NaT datetime columns default to an unzoned TimeZone (''), which
% MATLAB refuses to vertcat against a zoned datetime column. Used to keep
% every empty-fallback table concatenation-compatible with the real
% (zoned) data it may later be vertcat-ed against, and to reconstruct the
% TimeZone that plain CSV text does not carry (see read_pipeline_csv.m).

    for name = T.Properties.VariableNames
        if isdatetime(T.(name{1}))
            T.(name{1}).TimeZone = tz;
        end
    end
end
