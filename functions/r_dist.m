% R_DIST  integrated path distance between times a and b
% 
%     distance = r_dist(a, b, Pos_x, Pos_y, f_sample)
% 
% computes path length travelled by the mouse between times a and b 
% integrates small step-by-step movements between consecutive samples
% 
% INPUT
% a: start time in s 
% b: end time in s
% Pos_x: x-position of the mouse at each sample
% Pos_y: y-position of the mouse at each sample
% f_sample: sampling frequency
%
% NOTE
%   this is a refactored version of the original method:
%       function distance = r_dist(app,a,b,Pos_X,Pos_y,f_sample)
%   used inside the actimetry part of the app

function distance = r_dist(a, b, Pos_x, Pos_y, f_sample)

    a = a(:);
    b = b(:);
    Pos_x = Pos_x(:);
    Pos_y = Pos_y(:);

    n_samples = numel(Pos_x);
    n_intervals = numel(a);

    % preallocate output
    distance = zeros(n_intervals, 1);

    for k = 1:n_intervals
        ak = a(k);
        bk = b(k);

        n_steps = floor(f_sample * (bk - ak));
        if n_steps <= 0
            distance(k) = 0;
            continue;
        end 

        idx0 = round(ak * f_sample);
        idx0 = max(1 min(n_samples - 1, idx 0));

        d = 0;
        for step = 0:(n_steps-1)
            i1 = idx0 + step;
            i2 = i1 + 1;

            if i2 > n_samples
                break;
            end 

            dx = Pos_x(i2) - Pos_x(i1);
            dy = Pos_y(i2) - Pos_y(i1);
            d = d + sqrt(dx.^2 + dy.^2);
        end 

        distance(k) = d;
    end
end


