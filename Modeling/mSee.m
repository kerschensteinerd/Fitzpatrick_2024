function [mono] = mSee(fft_fshift,otfx,xl,yl)
% mSee - Apply optical transfer function to image in frequency domain
%
% This function applies a 1D OTF to a 2D image by:
% 1. Interpolating the 1D OTF to create a radially symmetric 2D OTF
% 2. Multiplying the image FFT by the 2D OTF
% 3. Converting back to spatial domain
% 4. Applying inverse gamma correction
%
% Inputs:
%   fft_fshift - 2D FFT of image (fftshift applied)
%   otfx       - 1D optical transfer function (at spatial frequencies xl)
%   xl         - Spatial frequency vector (cycles per degree) for OTF
%   yl         - Spatial frequency vector for second dimension
%
% Output:
%   mono       - Filtered image with gamma correction applied (0-1 range)
%
% Author: Fitzpatrick et al., 2024
    
    %create coordinates for 2d
    [Row,Col] = meshgrid(yl,xl);
    RCcoords(:,:,1) = Row;
    RCcoords(:,:,2) = Col;

    % Calculate distance from origin in frequency space
    dist = sqrt((RCcoords(:, :, 1) .^ 2 + RCcoords(:, :, 2) .^ 2));
    maxDist = min(max(yl),max(xl));
    dist_idx = (dist > maxDist);

    % Interpolate 1D OTF to 2D (radially symmetric)
    otfq = interp1(xl, otfx, dist, 'spline');
    otfq(dist_idx) = 0;  % Zero out beyond valid range
    % Create full 2D OTF (4 quadrants)
    otf2d = [flip([flip(otfq,2) otfq],1);flip(otfq,2) otfq];

    % Apply OTF in frequency domain
    otf_fftshift = fft_fshift.*otf2d;

    % Convert back to spatial domain
    otf_fft = ifftshift(otf_fftshift);
    lin = abs(ifft2(otf_fft));

    % Apply inverse gamma correction (sRGB standard)
    mono = zeros(size(lin));
    for i = 1:numel(lin)
        if lin(i)<=0.0031308
            mono(i) = lin(i)*12.92;  % Linear portion
        else
            mono(i) = (1.055*(lin(i))^(1/2.4)) - 0.055;  % Gamma portion
        end
    end
    
end