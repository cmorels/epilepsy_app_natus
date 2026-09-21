% DIST  to compute distance between mouse position at times a and b 
% 
%   distance = dist(a, b, Pos_x, Pos_y, f_sample)
% 
% INPUT
%   a: time in s
%   b: time in s (same as a)
%   Pos_x: x-position of the mouse at each sample
%   Pos_y: y-position of the mouse at each sample
%   f_sample: sampling frequency of position
% 
% NOTE
%   this is a refactored version of the original method:
%       function distance = dist(app,a,b,Pos_X,Pos_y,f_sample)
%   used inside the actimetry of the app

function distance = dist(a, b, Pos_x, Pos_y, f_sample)

    a = a(:);
    b = b(:);
    
    Pos_x = Pos_x(:);
    Pos_y = Pos_y(:);
    
    % Calculate the indices corresponding to the times a and b
    idx_a = round(a * f_sample);
    idx_b = round(b * f_sample);
    
    n = numel(Pos_x)
    idx_a = max(1, min(n, idx_a));
    idx_b = max(1, min(n, idx_b));
    
    % compute distance (euclidean)
    dx = Pos_x(idx_b) - Pos_x(idx_a);
    dy = Pos_y(idx_b) - Pos_y(idx_a);
    distance = sqrt(dx^2 + dy^2);

end