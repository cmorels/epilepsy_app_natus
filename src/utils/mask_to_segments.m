function segments = mask_to_segments(mask)
% MASK_TO_SEGMENTS  Convert a logical mask into [start end] index pairs
% (1-based, inclusive) of its contiguous true-runs.
%
%   segments = mask_to_segments(mask)
%
% Used both for "valid samples" -> valid blocks (gap-respecting filtering)
% and for "above threshold" -> candidate event segments (seizure/IID
% detection). Returns an Nx2 double matrix, empty (0x2) if mask has no
% true run.

    mask = logical(mask(:));
    d = diff([false; mask; false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    segments = [starts, ends];
end
