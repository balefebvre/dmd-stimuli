function [ A ] = moving_bars_generate_image( ws, hs, xb, yb, ab, wb, lb, ps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     ws: screen width (in px)
%     hs: screen height (in px)
%     x: bar x-coordinate (in µm)
%     y: bar y-coordinate (in µm)
%     a: bar angle (in rad)
%     wb: bar width (in µm)
%     lb: bar length (in µm)
%     ps: pixel_size (in µm)
    
    % Compute the position of each pixel.
    x = ((1:ws)' - (ws + 1) / 2) * ps;
    y = ((1:hs)' - (hs + 1) / 2) * ps;
    [X, Y] = meshgrid(x, y);
    
    % Find the pixels which compose the bar.
    nb_pixels = ws * hs;
    p = [X(:), Y(:)];
    p = p - repmat([xb, yb], nb_pixels, 1);
    a1 = p * [+cos(ab); +sin(ab)];
    a2 = p * [-sin(ab); +cos(ab)];
    b1 = abs(a1) < (0.5 * lb);
    b2 = abs(a2) < (0.5 * wb);
    b = b1 & b2;
    
    % Generate the image.
    A = zeros(hs, ws, 'uint8');
    A(b) = 1.0;
    A = flipud(A);
    
    return
    
end

