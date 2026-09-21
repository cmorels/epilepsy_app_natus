function y = filter_by_blocks(x, blocks, filter_fn)
% FILTER_BY_BLOCKS  Apply filter_fn independently to each contiguous block
% of x, never letting it see across a block boundary (i.e. never across a
% NaN gap). Positions outside every block are NaN in the output.
%
%   y = filter_by_blocks(x, blocks, filter_fn)
%
%   x         : numeric column vector (may contain NaN outside blocks)
%   blocks    : Nx2 [start end] index pairs (1-based, inclusive), e.g.
%               from mask_to_segments(~isnan(x))
%   filter_fn : function handle, y_block = filter_fn(x_block), applied to
%               each x(start:end) slice; must return a vector the same
%               length as its input.

    x = x(:);
    y = NaN(size(x));

    for i = 1:size(blocks, 1)
        s = blocks(i, 1);
        e = blocks(i, 2);
        y_block = filter_fn(x(s:e));
        if numel(y_block) ~= (e - s + 1)
            error('filter_by_blocks:BadFilterOutput', ...
                'filter_fn must return a vector the same length as its input block (block %d: in=%d, out=%d).', ...
                i, e - s + 1, numel(y_block));
        end
        y(s:e) = y_block(:);
    end
end
