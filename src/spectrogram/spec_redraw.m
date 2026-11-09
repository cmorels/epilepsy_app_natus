function files = spec_redraw(target, varargin)
% SPEC_REDRAW  Regenerate .fig / .png of finished events from their .mat
% only (no txt read, nothing recomputed) -- after a figure change, or when
% a PNG could not be written.
%
%   spec_redraw('...\06_spectrograms\073-s_20260706_130048\sz002_both.mat')
%   spec_redraw('...\073-s\06_spectrograms')          % every event .mat below it
%   spec_redraw(target, 'clim', [-10 10], 'freq_scale', 'log', ...)   % display overrides
%
% Name/value pairs override res.cfg display fields (clim, cmap,
% freq_scale, png_dpi, save_fig, save_png) for this redraw only.

    if isfolder(target)
        d = dir(fullfile(target, '**', '*.mat'));
        d = d(~contains({d.folder}, '_cache'));
        files = fullfile({d.folder}, {d.name});
    else
        files = {char(target)};
    end
    for i = 1:numel(files)
        S = load(files{i}, 'res');
        if ~isfield(S, 'res')
            continue;
        end
        res = S.res;
        for k = 1:2:numel(varargin)
            res.cfg.(varargin{k}) = varargin{k + 1};
        end
        res.cfg.visible = false;
        if strcmp(res.kind, 'seizure')
            fig = spec_plot_seizure(res);
        else
            fig = spec_plot_iid(res);
        end
        base = files{i}(1:end - 4);
        if res.cfg.save_fig
            set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
            savefig(fig, [base '.fig']);
        end
        if res.cfg.save_png
            spec_write_outputs('png', fig, [base '.png'], res.cfg.png_dpi);
        end
        close(fig);
        fprintf('redrawn %s\n', base);
    end
end
