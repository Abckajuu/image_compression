function show_gray(M, lims)
%SHOW_GRAY  Display a grayscale matrix with imshow() (imagesc as a fallback).
%   SHOW_GRAY(M) shows values 0..255. SHOW_GRAY(M, [lo hi]) uses custom limits.

    if nargin < 2, lims = [0 255]; end
    if exist('imshow', 'file')
        imshow(M, lims);
    else
        imagesc(M, lims); colormap(gray); axis image off;
    end
end
