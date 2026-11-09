function k = rp_first_sample(start, t, fs)
% RP_FIRST_SAMPLE  1-based index of the first sample at or after time t of
% a continuous recording whose sample 1 is at start (integer ms arithmetic,
% so a cut never moves by a rounding error).
    off_ms = int64(round(milliseconds(t - start)));
    k = double(idivide(off_ms * int64(fs) + 999, int64(1000), 'floor')) + 1;
end
