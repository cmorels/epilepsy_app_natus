function st = ll_status_style(status)
% LL_STATUS_STYLE  Shading style per ll_status, shared by every seizure
% figure (panorama, individual and joint event figures) so one event looks
% the same in every view.
%   'accepted'      red, solid edge
%   'in_band'       orange, dashed edge (review band)
%   'no_detection'  grey-blue, dotted edge (this channel did not detect it)
    switch status
        case 'accepted'
            st = struct('face', [1 0.7 0.7], 'edge', [0.85 0 0], 'line', '-', 'label', 'accepted');
        case 'in_band'
            st = struct('face', [1 0.85 0.55], 'edge', [0.9 0.45 0], 'line', '--', 'label', 'in\_band (review)');
        otherwise
            st = struct('face', [0.8 0.85 0.95], 'edge', [0.35 0.45 0.7], 'line', ':', 'label', 'no detection in this channel');
    end
end
