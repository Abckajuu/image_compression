%% MAIN  Image compression using Singular Value Decomposition (SVD)
%  Run this file. It loads an image, factorises it with svd(), rebuilds it
%  with different numbers of singular values, and reports storage and error.
clear; clc; close all;

%% 1. Settings -----------------------------------------------------------
imageFile  = 'peppers.png';     % any image on the MATLAB path, or your own file
fallback   = 'images/sample.png';
ks         = [5 10 20 50 100];  % ranks (numbers of singular values) to test
psnrTarget = 30;                % dB, a common rule of thumb for "acceptable" quality
here       = fileparts(mfilename('fullpath'));
outDir     = fullfile(here, 'results');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% 2. Load the image and represent it as a matrix ------------------------
try
    A = load_gray_image(imageFile);
    source = imageFile;
catch
    source = fullfile(here, fallback);
    fprintf('Could not read "%s". Using %s instead.\n', imageFile, fallback);
    A = load_gray_image(source);
end
[m, n] = size(A);
ks = ks(ks <= min(m, n));
fprintf('Image: %s\nSize : %d x %d  (%d numbers)\n', source, m, n, m * n);
fprintf('Rank : %d (full rank would be %d)\n\n', rank(A), min(m, n));

%% 3. Singular value decomposition: A = U * S * V' -----------------------
[U, S, V] = svd(A, 'econ');
s = diag(S);                                  % singular values, largest first
fprintf('Largest singular value %.1f, smallest %.3g\n', s(1), s(end));

% Link to eigenvalues: sigma_i^2 are the eigenvalues of A'*A
if min(m, n) <= 1000
    if m <= n, G = A * A'; else, G = A' * A; end
    ev = sort(eig((G + G') / 2), 'descend');
    nc = min(5, numel(s));
    gap = max(abs(sqrt(max(ev(1:nc), 0)) - s(1:nc)));
    fprintf('Top %d singular values vs sqrt(eigenvalues of A''A): max difference %.2e\n\n', nc, gap);
end

%% 4. Error for every rank k, straight from the singular values ----------
% Eckart-Young: ||A - A_k||_F^2 = sum of the squared discarded singular values
kAll      = (1:numel(s))';
total     = sum(s .^ 2);
retained  = cumsum(s .^ 2);
energyPct = 100 * retained / total;
relErrAll = sqrt(max(total - retained, 0) / total);
mseAll    = max(total - retained, 0) / (m * n);
psnrAll   = 10 * log10(255^2 ./ mseAll);
psnrAll(~isfinite(psnrAll)) = NaN;
storageAll = kAll * (m + n + 1);

%% 5. Reconstruct with k singular values and measure ---------------------
nk = numel(ks);
recon = cell(1, nk);
res   = cell(1, nk);
for i = 1:nk
    Ak = rank_k_approx(U, s, V, ks(i));
    Ak = min(max(Ak, 0), 255);                % keep valid pixel range
    recon{i} = Ak;
    res{i}   = compression_metrics(A, Ak, ks(i));
end

fprintf('%5s %10s %8s %8s %9s %9s %8s %8s\n', ...
        'k', 'Storage', 'Ratio', 'Saved%', 'RelErr', 'Theory', 'PSNR', 'Energy%');
for i = 1:nk
    r = res{i};
    fprintf('%5d %10d %7.2f:1 %7.1f%% %9.4f %9.4f %6.1fdB %7.2f%%\n', ...
            r.k, r.compressedStorage, r.ratio, r.savingsPct, r.relError, ...
            relErrAll(r.k), r.psnr, energyPct(r.k));
end
fprintf('(original storage: %d numbers; "Theory" ignores clipping to 0..255)\n\n', m * n);

%% 6. Summary: how many singular values are enough? ----------------------
breakEven = (m * n) / (m + n + 1);
fprintf('Compression only saves space when k < %.1f\n', breakEven);
for p = [90 95 99]
    k = find(energyPct >= p, 1);
    fprintf('%d%% of the energy needs k = %3d  (%.1f:1)\n', p, k, (m * n) / storageAll(k));
end
fprintf('(Energy is dominated by the average brightness, so judge quality by PSNR and by eye.)\n');
kq = find(psnrAll >= psnrTarget, 1);
if isempty(kq)
    fprintf('PSNR of %d dB is not reached with fewer than full rank.\n', psnrTarget);
else
    fprintf('PSNR >= %d dB first reached at k = %d  (%.1f:1, %.1f%% saved)\n', ...
            psnrTarget, kq, (m * n) / storageAll(kq), 100 * (1 - storageAll(kq) / (m * n)));
end

%% 7. Figures ------------------------------------------------------------
% Figure 1: original, compressed images, and where the error is
f1 = figure('Name', 'Reconstructions', 'Position', [40 40 min(1500, 250 * (nk + 1)) 560]);
subplot(2, nk + 1, 1); show_gray(A); title('Original');
subplot(2, nk + 1, nk + 2); axis off;
text(0.5, 0.5, sprintf('Top row: images\nBottom row:\nabsolute error'), 'HorizontalAlignment', 'center');
for i = 1:nk
    subplot(2, nk + 1, i + 1); show_gray(recon{i});
    title(sprintf('k = %d\n%.1f:1, %.1f dB', ks(i), res{i}.ratio, res{i}.psnr));
    subplot(2, nk + 1, nk + 2 + i); show_gray(abs(A - recon{i}), [0 64]);
    title(sprintf('error, k = %d', ks(i)));
end

% Figure 2: singular value spectrum and cumulative energy
f2 = figure('Name', 'Singular values', 'Position', [60 60 1000 400]);
subplot(1, 2, 1);
semilogy(kAll, s, 'LineWidth', 1.3); grid on;
xlabel('Index i'); ylabel('Singular value \sigma_i'); title('Singular value spectrum');
subplot(1, 2, 2);
plot(kAll, energyPct, 'LineWidth', 1.3); hold on; grid on;
k95 = find(energyPct >= 95, 1);
plot(k95, energyPct(k95), 'ro', 'MarkerFaceColor', 'r');
text(k95, energyPct(k95) - 6, sprintf('  95%% at k = %d', k95));
xlabel('Number of singular values k'); ylabel('Energy retained (%)');
title('Cumulative energy'); ylim([0 100]);

% Figure 3: error analysis and storage
f3 = figure('Name', 'Error analysis', 'Position', [80 80 1300 400]);
subplot(1, 3, 1);
semilogy(kAll, relErrAll, 'LineWidth', 1.3); hold on; grid on;
semilogy(ks, cellfun(@(r) r.relError, res), 'ro', 'MarkerFaceColor', 'r');
xlabel('k'); ylabel('Relative error ||A - A_k||_F / ||A||_F'); title('Reconstruction error');
subplot(1, 3, 2);
plot(kAll, psnrAll, 'LineWidth', 1.3); hold on; grid on;
plot([1 numel(s)], [psnrTarget psnrTarget], 'k--');
plot(ks, cellfun(@(r) r.psnr, res), 'ro', 'MarkerFaceColor', 'r');
xlabel('k'); ylabel('PSNR (dB)'); title(sprintf('Quality (dashed: %d dB)', psnrTarget));
subplot(1, 3, 3);
plot(kAll, storageAll, 'LineWidth', 1.3); hold on; grid on;
plot([1 numel(s)], [m * n m * n], 'r--');
xlabel('k'); ylabel('Numbers stored'); title('Storage: compressed vs original (dashed)');

saveas(f1, fullfile(outDir, '1_reconstructions.png'));
saveas(f2, fullfile(outDir, '2_singular_values.png'));
saveas(f3, fullfile(outDir, '3_error_analysis.png'));
fprintf('\nFigures saved in %s\n', outDir);
