function [sz, man, ids] = cs_iid_channel_zones(zones, ch, buf)
% CS_IID_CHANNEL_ZONES  The zones of cs_iid_zones.m that apply to one
% channel (same mouse, overlapping the channel's span once buffered),
% converted to the file's own time base (t_rel = POSIX - session start, s).
%
%   zones : cs_iid_zones.m's Z.zones
%   ch    : one row of ch_tbl (mouse, start_posix, duration_s)
%   buf   : cfg.iid_excl_buffer_s
%   sz    : table (id = zone_id, start_s, end_s) WITHOUT buffer -> detect_iid's
%           seizures argument (detect_iid adds cfg.iid.exclusion_buffer_s)
%   man   : Nx2 [start_s end_s] WITH +/- buf -> cfg.iid.exclusion_zones_manual
%   ids   : zone_id of every zone used (seizure zones first, then manual)

    t0 = ch.start_posix;
    dur = ch.duration_s;
    k = strcmp(zones.mouse, ch.mouse{1}) & zones.end_posix + buf > t0 & zones.start_posix - buf < t0 + dur;
    z = zones(k, :);
    s = reshape(z.start_posix - t0, [], 1);
    e = reshape(z.end_posix - t0, [], 1);
    id = reshape(z.zone_id, [], 1);
    isz = reshape(strcmp(z.pass_as, 'seizures'), [], 1);
    col = @(v, k) reshape(v(k), [], 1);   % a scalar indexed by a false scalar is 0x0, not 0x1
    sz = table(col(id, isz), col(s, isz), col(e, isz), 'VariableNames', {'id', 'start_s', 'end_s'});
    man = [col(s, ~isz) - buf, col(e, ~isz) + buf];
    ids = [col(id, isz); col(id, ~isz)];
end
