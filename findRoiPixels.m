function [pixels,pixelIndices] = findRoiPixels(center,radius,imageWidth,imageHeight)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com
%
% Return Nx2 array of [x, y] pixels in ROI.

%% Init.
pixels = [];
pixelIndices = [];

%% Find pixels in ROI.
cx = center(1);
cy = center(2);
x0 = round(cx);
y0 = round(cy);
pixels = [x0, y0];
r = ceil(radius);
for x = x0-r:x0+r
    for y = y0-r:y0+r
        if (x ~= x0) || (y ~= y0)
            dx = x-x0;
            dy = y-y0;
            dr = sqrt(dx^2 + dy^2);
            if dr <= radius
                pixels = [pixels;[x,y]];
            end
        end
    end
end
pixels = uint32(pixels);

%% Remove pixels outside of image boundaries.
x = pixels(:,1);
y = pixels(:,2);
outX = union(find(x < 1),find(x > imageWidth));
outY = union(find(y < 1),find(y > imageHeight));
outsideImagePixels = union(outX,outY);
if ~isempty(outsideImagePixels)
    pixels(outsideImagePixels,:) = [];
end

%% Convert pixels to pixel matrix indices.
numPixels = size(pixels,1);
if numPixels > 0
    pixelIndices = zeros(1,numPixels);
    for i = 1:numPixels
        x = pixels(i,1);
        y = pixels(i,2);
        row = y;
        col = x;
        index = (col-1)*imageHeight+row;
        pixelIndices(i) = index;
    end
    pixelIndices = uint32(pixelIndices);
end

end
