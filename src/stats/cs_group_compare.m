function g = cs_group_compare(e, c, cfg)
% CS_GROUP_COMPARE  Epileptic (e) vs control (c) per-mouse values: fold
% change of medians with a stratified bootstrap percentile CI
% (cfg.n_boot, cfg.rng_seed), Mann-Whitney p and Cliff's delta
% (P(e > c) - P(e < c)). Same computation as cs_iid_stats.m's group table.

    g.fold = median(e) / median(c);
    rng(cfg.rng_seed);
    fb = zeros(cfg.n_boot, 1);
    ne = numel(e); nc = numel(c);
    for b = 1:cfg.n_boot
        fb(b) = median(e(randi(ne, ne, 1))) / median(c(randi(nc, nc, 1)));
    end
    fb = fb(~isnan(fb));
    g.ci = prctile(fb, [2.5 97.5]);
    g.p = ranksum(e, c);
    d = sign(e(:) - c(:)');
    g.cliff = mean(d(:));
end
