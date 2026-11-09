function spec_log(varargin)
% SPEC_LOG  Print a line and append it to the run log.
%
%   spec_log('open', log_path)   start appending to log_path
%   spec_log('close')
%   spec_log(fmt, args...)       fprintf-style line (timestamped in the file)

    persistent fid
    if nargin >= 1 && ischar(varargin{1})
        switch varargin{1}
            case 'open'
                if ~isempty(fid) && fid > 0
                    fclose(fid);
                end
                fid = fopen(varargin{2}, 'a', 'n', 'UTF-8');
                return;
            case 'close'
                if ~isempty(fid) && fid > 0
                    fclose(fid);
                end
                fid = [];
                return;
        end
    end
    msg = sprintf(varargin{:});
    fprintf('%s\n', msg);
    if ~isempty(fid) && fid > 0
        fprintf(fid, '%s  %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), msg);
    end
end
