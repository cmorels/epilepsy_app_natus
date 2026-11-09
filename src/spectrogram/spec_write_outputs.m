function out = spec_write_outputs(action, varargin)
% SPEC_WRITE_OUTPUTS  Files of the spectrogram module (stage 06). Nothing
% is ever written in 01_txt ... 05_summaries.
%
%   paths = spec_write_outputs('event', res, fig, out_dir, git)
%       <out_dir>\<id>_<region>.mat  : res (everything needed to redraw:
%                                      spec_plot_seizure / spec_plot_iid),
%                                      features, kind, git
%       <out_dir>\<id>_<region>.fig / .png (png at cfg.spec.png_dpi)
%   T = spec_write_outputs('summary', output_root, summary_dir)
%       Rebuilds summary_dir\spec_features_seizures.csv and
%       spec_features_iid.csv from EVERY event .mat under
%       <output_root>\<subject>\06_spectrograms\ (so re-running one event
%       replaces its row and earlier events are kept).

    switch action
        case 'event'
            out = write_event(varargin{:});
        case 'summary'
            out = write_summary(varargin{:});
        otherwise
            error('spec_write_outputs:BadAction', 'Unknown action ''%s''.', action);
    end
end

%% ======================================================================
function paths = write_event(res, fig, out_dir, git)
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    base = fullfile(out_dir, sprintf('%s_%s', safe(res.id), res.region_tag));
    features = res.features;
    features.git_head = repmat({git.head}, height(features), 1);
    features.git_dirty = repmat(git.dirty, height(features), 1);
    features.processed_at = repmat({char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'))}, height(features), 1);
    paths = struct('mat', [base '.mat'], 'fig', '', 'png', '');
    features.mat_file = repmat({paths.mat}, height(features), 1);
    features.png_file = repmat({''}, height(features), 1);
    if res.cfg.save_png
        paths.png = [base '.png'];
        features.png_file(:) = {paths.png};
    end
    res.features = features;
    kind = res.kind;
    save(paths.mat, 'res', 'features', 'kind', 'git', '-v7.3');
    if ~isempty(fig) && isgraphics(fig)
        if res.cfg.save_fig
            paths.fig = [base '.fig'];
            set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');   % hidden while batch-saving, visible when opened
            savefig(fig, paths.fig);
        end
        if res.cfg.save_png
            exportgraphics(fig, paths.png, 'Resolution', res.cfg.png_dpi);
        end
    end
end

function T = write_summary(output_root, summary_dir)
    if ~isfolder(summary_dir)
        mkdir(summary_dir);
    end
    files = dir(fullfile(output_root, '*', '06_spectrograms', '**', '*.mat'));
    files = files(~contains({files.folder}, '_cache'));
    parts = struct('seizure', {{}}, 'iid', {{}});
    for i = 1:numel(files)
        p = fullfile(files(i).folder, files(i).name);
        try
            S = load(p, 'features', 'kind');
        catch
            continue;
        end
        if isfield(S, 'kind') && isfield(S, 'features') && isfield(parts, S.kind)
            parts.(S.kind){end+1} = S.features;
        end
    end
    T = struct();
    for k = {'seizure', 'iid'}
        f = fullfile(summary_dir, sprintf('spec_features_%s.csv', plural(k{1})));
        if isempty(parts.(k{1}))
            T.(k{1}) = table();
            continue;
        end
        T.(k{1}) = stack(parts.(k{1}));
        writetable(T.(k{1}), f);
    end
end

function s = plural(k)
    if strcmp(k, 'seizure')
        s = 'seizures';
    else
        s = 'iid';
    end
end

function T = stack(parts)
% vertcat tables whose columns may differ (older outputs): union of columns.
    names = {};
    for i = 1:numel(parts)
        names = [names, setdiff(parts{i}.Properties.VariableNames, names, 'stable')]; %#ok<AGROW>
    end
    for i = 1:numel(parts)
        miss = setdiff(names, parts{i}.Properties.VariableNames, 'stable');
        for m = miss
            parts{i}.(m{1}) = repmat({''}, height(parts{i}), 1);
        end
        parts{i} = parts{i}(:, names);
    end
    T = vertcat(parts{:});
end

function s = safe(s)
    s = regexprep(char(s), '[^\w\-]', '_');
end
