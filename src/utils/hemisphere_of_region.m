function h = hemisphere_of_region(region, cfg)
% HEMISPHERE_OF_REGION  'left' | 'right' | 'unknown', from the REGION name
% only (cfg.output.hemisphere.left_regions / right_regions, case-insensitive,
% trimmed). Never from the channel label: physical channels change between
% recordings and animals (A5C2/A5C4, A7C1/A7C3, ...), the region name set by
% the recording log does not.
%
%   h = hemisphere_of_region(region, cfg)     region: char or cellstr
    if iscell(region)
        h = cellfun(@(r) hemisphere_of_region(r, cfg), region, 'UniformOutput', false);
        return;
    end
    r = lower(strtrim(char(region)));
    if any(strcmp(r, lower(strtrim(cfg.output.hemisphere.left_regions))))
        h = 'left';
    elseif any(strcmp(r, lower(strtrim(cfg.output.hemisphere.right_regions))))
        h = 'right';
    else
        h = 'unknown';
    end
end
