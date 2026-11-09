function E = rp_edf_read(path, labels)
% RP_EDF_READ  Header and DIGITAL samples (int16, untouched) of an EDF/EDF+.
%
%   E = rp_edf_read(path)          header only
%   E = rp_edf_read(path, labels)  plus the samples of these signals
%
% E.hdr   : version, patient, recording, startdate, starttime, header_bytes,
%           reserved, nrec, rec_dur_s, ns
% E.sig   : struct array, one per signal: label, transducer, dim, pmin, pmax,
%           dmin, dmax, prefilter, spr, reserved, fs
% E.x     : nsamp x numel(labels) int16 (labels order), concatenated records
% E.ann   : (2*spr) x nrec uint8, raw bytes of 'EDF Annotations' (if any)
% E.tal0_s: onset of the first record from its time-keeping TAL (0 if none)
% phys = pmin + (double(x) - dmin) * (pmax - pmin) / (dmax - dmin)  (rp_edf_phys)

    fid = fopen(path, 'r', 'ieee-le');
    if fid < 0, error('rp_edf_read:Open', 'Cannot open %s', path); end
    cl = onCleanup(@() fclose(fid));
    h = fread(fid, 256, '*char')';
    E.hdr = struct('version', strtrim(h(1:8)), 'patient', strtrim(h(9:88)), 'recording', strtrim(h(89:168)), ...
        'startdate', h(169:176), 'starttime', h(177:184), 'header_bytes', str2double(h(185:192)), ...
        'reserved', strtrim(h(193:236)), 'nrec', str2double(h(237:244)), 'rec_dur_s', str2double(h(245:252)), ...
        'ns', str2double(h(253:256)));
    ns = E.hdr.ns;
    H = fread(fid, ns * 256, '*char')';
    f = @(off, w) arrayfun(@(i) strtrim(H(off * ns + (i - 1) * w + (1:w))), 1:ns, 'UniformOutput', false);
    lab = f(0, 16); tr = f(16, 80); dim = f(96, 8); pmin = f(104, 8); pmax = f(112, 8);
    dmin = f(120, 8); dmax = f(128, 8); pf = f(136, 80); spr = f(216, 8); res = f(224, 32);
    E.sig = struct('label', lab, 'transducer', tr, 'dim', dim, 'pmin', num2cell(str2double(pmin)), ...
        'pmax', num2cell(str2double(pmax)), 'dmin', num2cell(str2double(dmin)), 'dmax', num2cell(str2double(dmax)), ...
        'prefilter', pf, 'spr', num2cell(str2double(spr)), 'reserved', res);
    for i = 1:ns, E.sig(i).fs = E.sig(i).spr / E.hdr.rec_dur_s; end
    d = dir(path);
    sz = d.bytes;
    rec_bytes = 2 * sum([E.sig.spr]);
    if sz ~= E.hdr.header_bytes + E.hdr.nrec * rec_bytes
        error('rp_edf_read:Size', '%s: %d bytes, header says %d', path, sz, E.hdr.header_bytes + E.hdr.nrec * rec_bytes);
    end
    off = [0 cumsum(2 * [E.sig.spr])];

    E.x = zeros(0, 0, 'int16');
    if nargin > 1
        E.x = zeros(E.hdr.nrec * E.sig(find(strcmp({E.sig.label}, labels{1}), 1)).spr, numel(labels), 'int16');
        for k = 1:numel(labels)
            i = find(strcmp({E.sig.label}, labels{k}), 1);
            if isempty(i), error('rp_edf_read:Label', '%s: no signal "%s"', path, labels{k}); end
            fseek(fid, E.hdr.header_bytes + off(i), 'bof');
            v = fread(fid, E.sig(i).spr * E.hdr.nrec, sprintf('%d*int16=>int16', E.sig(i).spr), rec_bytes - 2 * E.sig(i).spr);
            E.x(:, k) = v;
        end
    end

    E.ann = zeros(0, 0, 'uint8');
    E.tal0_s = 0;
    ia = find(strcmp({E.sig.label}, 'EDF Annotations'), 1);
    if ~isempty(ia)
        fseek(fid, E.hdr.header_bytes + off(ia), 'bof');
        n = 2 * E.sig(ia).spr;
        E.ann = reshape(fread(fid, n * E.hdr.nrec, sprintf('%d*uint8=>uint8', n), rec_bytes - n), n, E.hdr.nrec);
        tok = regexp(char(E.ann(:, 1)'), '^([+-][0-9.]+)\x14\x14', 'tokens', 'once');
        if ~isempty(tok), E.tal0_s = str2double(tok{1}); end
    end
end
