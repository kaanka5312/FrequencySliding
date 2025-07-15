%% --- helper that works on ONE matrix ----------------------------------
function out = row_acw(mat, fs)
    % mat : [nRows  x  nTimePoints]
    % fs  : sampling frequency (Hz)
    nRows        = size(mat,1);
    acw0         = zeros(nRows,1);
    acw50        = zeros(nRows,1);

    for r = 1:nRows                     % one ACW per row
        [acw0(r), acw50(r)] = acw(mat(r,:), fs, false);
    end

    % pack any way you like ─ here a 2-column matrix
    out = [acw0 acw50];                 % [nRows × 2]
end