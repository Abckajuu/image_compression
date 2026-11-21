# Image Compression Using Singular Value Decomposition (SVD)

A MATLAB mini project that compresses a grayscale image by keeping only the largest singular values of its matrix, then measures how much storage is saved and how much quality is lost.

**MATLAB concepts:** matrix operations, rank, eigenvalues, singular values, `svd()`, error analysis.
**Syllabus mapping:** Units 1, 2 and 5.

## How it works

A grayscale image is an m x n matrix A. Its SVD is `A = U * S * V'`, where S holds the singular values in decreasing order. Keeping only the first k terms gives

```
A_k = sum_{i=1..k} sigma_i * u_i * v_i'
```

By the Eckart-Young theorem, A_k is the best rank-k approximation of A. Instead of m*n numbers you store `k*(m+n+1)` numbers (k columns of U, k of V, and k singular values), so compression helps only while `k < m*n / (m+n+1)`.

The error is known without reconstructing anything: `||A - A_k||_F^2` equals the sum of the squared discarded singular values. The project checks this against the measured error. It also checks the link to eigenvalues: each `sigma_i^2` is an eigenvalue of `A'*A`.

## Run it

1. Open this folder in MATLAB.
2. Run `main`.

It reads `peppers.png` (shipped with MATLAB). If that file is not found it uses `images/sample.png`. To use your own picture, change `imageFile` at the top of `main.m`. You can also change the list of ranks `ks` and the quality target `psnrTarget`.

## Files

```
main.m                 runs the whole project and saves the figures
load_gray_image.m      imread + conversion to a grayscale matrix
rank_k_approx.m        builds A_k from U, the singular values and V
compression_metrics.m  storage, compression ratio, relative error, MSE, PSNR
show_gray.m            displays a matrix with imshow
images/sample.png      test image
results/               figures are saved here when you run main
```

## Output

The command window prints the image size, rank, a table for each k (storage, ratio, space saved, relative error, theoretical error, PSNR, energy kept) and a summary of how many singular values reach a given quality. Three figures are shown and saved to `results/`:

1. Original and compressed images for each k, with a row showing where the error is.
2. Singular value spectrum and cumulative energy.
3. Reconstruction error, PSNR and storage as k changes.

## Reading the results

- **PSNR** is the main quality measure. Around 30 dB looks acceptable to most eyes; above 40 dB is very hard to tell from the original.
- The **energy** percentage is dominated by the image's average brightness, so it reaches 90% with very few values. Judge quality by PSNR and by eye, not by energy alone.
- Smooth images compress with small k. Images with fine texture or noise need larger k.

## Testing note

The calculations were tested in GNU Octave on `images/sample.png`: the measured error matches the theoretical error from the singular values, and the top singular values match the square roots of the eigenvalues of `A'*A`. The plotting code could not be run in that test setup, so check the three figures when you first run `main` in MATLAB.
