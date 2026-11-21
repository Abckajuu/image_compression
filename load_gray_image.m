function G = load_gray_image(file)
%LOAD_GRAY_IMAGE  Read an image with imread() and return a grayscale matrix.
%   G = LOAD_GRAY_IMAGE(file) returns a double matrix with values in 0..255.
%   Works for RGB, grayscale and indexed images. It uses the standard
%   luminance weights, so the Image Processing Toolbox (rgb2gray) is not needed.

    [I, map] = imread(file);

    if ~isempty(map)                       % indexed image -> RGB
        I = ind2rgb(I, map) * 255;
    elseif isa(I, 'uint16')                % 16-bit image -> 8-bit range
        I = double(I) / 257;
    end
    I = double(I);

    if ndims(I) == 3 && size(I, 3) >= 3    % RGB -> gray
        G = 0.2989 * I(:,:,1) + 0.5870 * I(:,:,2) + 0.1140 * I(:,:,3);
    else
        G = I(:,:,1);
    end
end
