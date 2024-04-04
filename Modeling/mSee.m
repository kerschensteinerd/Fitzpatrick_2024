function [mono] = mSee(fft_fshift,otfx,xl,yl)
    
    %create coordinates for 2d
    [Row,Col] = meshgrid(yl,xl);
    RCcoords(:,:,1) = Row;
    RCcoords(:,:,2) = Col;

    dist = sqrt((RCcoords(:, :, 1) .^ 2 + RCcoords(:, :, 2) .^ 2));
    maxDist = min(max(yl),max(xl));
    dist_idx = (dist > maxDist);

    otfq = interp1(xl, otfx, dist, 'spline');
    otfq(dist_idx) = 0;
    otf2d = [flip([flip(otfq,2) otfq],1);flip(otfq,2) otfq];

    otf_fftshift = fft_fshift.*otf2d;

    otf_fft = ifftshift(otf_fftshift);

    lin = abs(ifft2(otf_fft));

    mono = zeros(size(lin));
    for i = 1:numel(lin)
        if lin(i)<=0.0031308
            mono(i) = lin(i)*12.92;
        else
            mono(i) = (1.055*(lin(i))^(1/2.4)) - 0.055;
        end
    end
    
end