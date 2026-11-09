function bad = cs_fig_dates(fig)
% CS_FIG_DATES  Every visible text of a figure (titles, subtitles, axis
% labels, tick labels, legends, text objects) that looks like a calendar
% date (yyyy-mm-dd, dd/mm, dd.mm.yy, a 20xx year or a month name). v2
% section 7 figures must have none. bad: cellstr of offending strings.

    s = {};
    for t = [findall(fig, 'Type', 'text'); findall(fig, 'Type', 'textboxshape')]'
        s = [s; cellstr(string(t.String(:)))]; %#ok<AGROW>
    end
    for ax = findall(fig, 'Type', 'axes')'
        s = [s; cellstr(string(ax.XTickLabel)); cellstr(string(ax.YTickLabel))]; %#ok<AGROW>
        for p = {'Title', 'Subtitle', 'XLabel', 'YLabel'}
            if isprop(ax, p{1}), s = [s; cellstr(string(ax.(p{1}).String))]; end %#ok<AGROW>
        end
    end
    for lg = findall(fig, 'Type', 'legend')'
        s = [s; cellstr(string(lg.String(:)))]; %#ok<AGROW>
    end
    pat = ['\d{4}-\d{1,2}-\d{1,2}|\d{1,2}/\d{1,2}|\d{1,2}\.\d{1,2}\.\d{2,4}|\<20\d{2}\>|' ...
        '\<(ene|feb|mar|abr|may|jun|jul|ago|sep|sept|oct|nov|dic|jan|apr|aug|dec)\>|' ...
        '\<(enero|febrero|marzo|abril|mayo|junio|julio|agosto|septiembre|octubre|noviembre|diciembre|' ...
        'january|february|march|april|june|july|august|september|october|november|december)\>'];
    hit = ~cellfun(@isempty, regexpi(s, pat, 'once'));
    bad = unique(s(hit));
end
