function rois = findRoisInMask(binaryMask,minRoiArea,minRoiSeparation)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com

% Find ROIs in image mask using regionprops.

% args:
% -----
% binaryMask:           Image mask in which to find ROIs using regionprops.
% minRoiArea:           Minimum ROI area.
% minRoiSeparation:     Minimum distance between ROI centroids.

% return:
% -------
% rois:                 Struct array of ROIs from regionprops after culling
%                       ROIs that are too small or too close together.

%% Init.
rois = [];

%% Make sure mask is binary.
binaryMask = (binaryMask > 0);

%% Find ROI properties.
rois = regionprops(binaryMask,'Centroid','Area','BoundingBox','PixelIdxList','PixelList');

%% Remove ROIs that are too small.
areas = vertcat(rois.Area);
badRois = find(areas < minRoiArea);
if ~isempty(badRois)
    rois(badRois) = [];
end

%% Remove ROIs that are too close together.
badRois = [];
for i = 1:length(rois)
    for j = 1:length(rois)
        if (i ~= j) && isempty(find(badRois == i)) && isempty(find(badRois == j))
            xi = rois(i).Centroid(1);
            yi = rois(i).Centroid(2);
            xj = rois(j).Centroid(1);
            yj = rois(j).Centroid(2);
            interRoiDistance = sqrt((xi-xj)^2 + (yi-yj)^2);
            if interRoiDistance < minRoiSeparation
                % Decide wether to throw out ROI i or j.
                if rois(i).Area >= rois(j).Area
                    badRois(end+1) = j;
                else
                    badRois(end+1) = i;
                end
            end
        end
    end
end
if ~isempty(badRois)
    rois(badRois) = [];
end

%% Plot mask with ROIs marked.
% figure;
% imshow(binaryMask); hold on;
% centroids = vertcat(rois.Centroid);
% plot(centroids(:,1),centroids(:,2),'r+');

end
