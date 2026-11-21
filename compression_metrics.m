function m = compression_metrics(A, Ak, k)
%COMPRESSION_METRICS  Storage and error figures for a rank-k approximation.
%   m = COMPRESSION_METRICS(A, Ak, k) returns a struct with:
%     originalStorage    numbers stored for A            (rows * cols)
%     compressedStorage  numbers stored for U_k, s_k, V_k   k*(rows + cols + 1)
%     ratio              originalStorage / compressedStorage
%     savingsPct         percentage of storage saved (negative = larger than A)
%     relError           ||A - Ak||_F / ||A||_F
%     mse, psnr          mean squared error and peak SNR in dB (peak = 255)

    [r, c] = size(A);
    m.k                 = k;
    m.originalStorage   = r * c;
    m.compressedStorage = k * (r + c + 1);
    m.ratio             = m.originalStorage / m.compressedStorage;
    m.savingsPct        = 100 * (1 - m.compressedStorage / m.originalStorage);

    E = A - Ak;
    m.relError = norm(E, 'fro') / norm(A, 'fro');
    m.mse      = mean(E(:) .^ 2);
    if m.mse == 0
        m.psnr = Inf;
    else
        m.psnr = 10 * log10(255^2 / m.mse);
    end
end
