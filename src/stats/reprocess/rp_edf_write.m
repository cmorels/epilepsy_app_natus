function rp_edf_write(path, W)
% RP_EDF_WRITE  Write an EDF+C file with DIGITAL samples copied as given
% (no rescaling) plus an 'EDF Annotations' signal.
%
% W.patient, W.recording : header text fields
% W.start                : datetime of the header start (whole second)
% W.frac_ns              : onset of record 1 after W.start, in ns (< 1e9);
%                          written as the time-keeping TAL of every record
% W.spr                  : samples per record (record duration = spr / fs = 1 s here)
% W.rec_dur_s            : record duration (s)
% W.sig                  : struct array (rp_edf_read's E.sig fields) of the data signals
% W.x                    : n x numel(W.sig) int16, n a multiple of W.spr
% W.tals                 : n_ann x 2 cell {onset_s (relative to W.start), text bytes uint8}
%                          ('\x14'-separated texts, as in the source TAL)
% W.ann_spr              : samples (2 bytes) reserved for annotations per record

    ns = numel(W.sig) + 1;
    n = size(W.x, 1);
    if mod(n, W.spr) ~= 0, error('rp_edf_write:Records', 'n = %d is not a multiple of spr = %d', n, W.spr); end
    nrec = n / W.spr;
    nb = 2 * W.ann_spr;

    % annotation bytes per record: time keeping TAL first, then event TALs
    A = zeros(nb, nrec, 'uint8');
    used = zeros(1, nrec);
    for r = 1:nrec
        b = uint8(sprintf('+%s%s%s%s', num2str(r - 1), frac_str(W.frac_ns), char([20 20 0])));
        A(1:numel(b), r) = b;
        used(r) = numel(b);
    end
    for a = 1:size(W.tals, 1)
        b = [uint8(['+' onset_str(W.tals{a, 1}) char(20)]), W.tals{a, 2}(:)', uint8(0)];
        if numel(b) > nb - 25, error('rp_edf_write:LongTAL', 'annotation of %d bytes does not fit', numel(b)); end
        r = min(max(floor(W.tals{a, 1} - W.frac_ns / 1e9) + 1, 1), nrec);
        while r <= nrec && used(r) + numel(b) > nb, r = r + 1; end
        if r > nrec, error('rp_edf_write:NoRoom', 'no room left for annotation %d', a); end
        A(used(r) + (1:numel(b)), r) = b;
        used(r) = used(r) + numel(b);
    end

    fid = fopen(path, 'w', 'ieee-le');
    if fid < 0, error('rp_edf_write:Open', 'Cannot create %s', path); end
    cl = onCleanup(@() fclose(fid));
    hdr = [pad_('0', 8), pad_(W.patient, 80), pad_(W.recording, 80), ...
        pad_(char(W.start, 'dd.MM.yy'), 8), pad_(char(W.start, 'HH.mm.ss'), 8), pad_(num2str(256 * (ns + 1)), 8), ...
        pad_('EDF+C', 44), pad_(num2str(nrec), 8), pad_(num2str(W.rec_dur_s), 8), pad_(num2str(ns), 4)];
    s = W.sig;
    lab = [{s.label}, {'EDF Annotations'}];
    fld = @(v, w) strjoin(cellfun(@(x) pad_(x, w), v, 'UniformOutput', false), '');
    num = @(v) cellfun(@(x) num2str(x), num2cell(v), 'UniformOutput', false);
    hdr = [hdr, fld(lab, 16), fld([{s.transducer}, {''}], 80), fld([{s.dim}, {''}], 8), ...
        fld([num([s.pmin]), {'-1'}], 8), fld([num([s.pmax]), {'1'}], 8), fld([num([s.dmin]), {'-32768'}], 8), ...
        fld([num([s.dmax]), {'32767'}], 8), fld([{s.prefilter}, {''}], 80), ...
        fld([num(repmat(W.spr, 1, numel(s))), {num2str(W.ann_spr)}], 8), fld(repmat({''}, 1, ns), 32)];
    if numel(hdr) ~= 256 * (ns + 1), error('rp_edf_write:Header', 'header has %d bytes', numel(hdr)); end
    fwrite(fid, hdr, 'char');

    D = zeros(numel(s) * W.spr + W.ann_spr, nrec, 'int16');
    for c = 1:numel(s)
        D((c - 1) * W.spr + (1:W.spr), :) = reshape(W.x(:, c), W.spr, nrec);
    end
    D(numel(s) * W.spr + 1:end, :) = reshape(typecast(A(:), 'int16'), W.ann_spr, nrec);
    fwrite(fid, D, 'int16');
end

function s = pad_(s, w)
    s = char(s);
    if numel(s) > w, error('rp_edf_write:Field', '"%s" longer than %d', s, w); end
    s = [s, repmat(' ', 1, w - numel(s))];
end

function s = frac_str(ns)
    if ns == 0, s = ''; return; end
    s = ['.' regexprep(sprintf('%09d', ns), '0+$', '')];
end

function s = onset_str(t)
    s = regexprep(sprintf('%.6f', t), '\.?0+$', '');
end
