function cs_write_sheet(file, sheet, T)
% CS_WRITE_SHEET  writetable to one sheet; datetimes become text (local
% wall clock, no zone ambiguity in Excel), POSIX helper columns are kept,
% an empty table writes a one-line note instead of nothing.

    if isempty(T) || height(T) == 0
        T = table({'(sin filas)'}, 'VariableNames', {'nota'});
    end
    for v = T.Properties.VariableNames
        x = T.(v{1});
        if isdatetime(x)
            if isempty(x.TimeZone) && all(timeofday(x(~isnat(x))) == 0)
                x.Format = 'yyyy-MM-dd';
            else
                x.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
            end
            T.(v{1}) = cellstr(string(x));
        elseif islogical(x)
            T.(v{1}) = double(x);
        end
    end
    writetable(T, file, 'Sheet', sheet, 'WriteMode', 'overwritesheet');
end
