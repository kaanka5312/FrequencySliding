%% ---------- helper (one matrix → one column vector) -------------------
function MF_col = row_mf(mat, fs, freqrange)
    nRows  = size(mat,1);
    MF_col = zeros(nRows,1);

    for r = 1:nRows
        MF_col(r) = mf(mat(r,:), fs, freqrange);
    end
end