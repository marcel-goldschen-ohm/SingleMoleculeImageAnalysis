function rois = findAndProjectRoisForZmwAlexFret(...
    fluorescenceImageStackFilePath,brightFieldImageFilePath,...
    minRoiArea,minRoiSeparation,acceptorRoiLocalSearch,donorRoiLocalSearch,roiRadius)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com
%
% Requires Coherent Point Drift (https://sites.google.com/site/myronenko/research/cpd) 
%   for cpd_register() function.
%
% args:
% -----
% fluorescenceImageStackFilePath:   File path of 1024x512xN TIFF image
%                                   stack with N frames where each frame
%                                   contains two 512x512 camera images
%                                   sie-by-side (left=acceptor,
%                                   right=donor).
% brightFieldImageFilePath:         File path of 1024x512 TIFF bright field
%                                   image of ZMW arrays from two camera
%                                   images placed side-by-side
%                                   (left=acceptor, right=donor).
%
% return:
% -------
% rois:                             Struct array of ROIs (e.g. as returned
%                                   by regionprops) for ZMWs containing
%                                   acceptor molecules.

%% Init.
if ~exist('minRoiArea','var')
    minRoiArea = 2;
end
if ~exist('minRoiSeparation','var')
    minRoiSeparation = 5;
end
if ~exist('acceptorRoiLocalSearch','var')
    acceptorRoiLocalSearch = 1;
end
if ~exist('donorRoiLocalSearch','var')
    donorRoiLocalSearch = 1;
end
if ~exist('roiRadius','var')
    roiRadius = 2.5;
end

%% Make sure we have an image stack to analyze.
if ~exist('fluorescenceImageStackFilePath','var')
    [file,path] = uigetfile('*.tif','Load dual camera fluorescence image stack.');
    if ~isequal(file,0)
        fluorescenceImageStackFilePath = [path file];
    else
        errordlg('ERROR: findAndProjectRoisForZmwAlexFret requires a valid TIFF image stack.');
        return
    end
end

%% Figure to visualize ROI detection.
roisFig = figure;

%% Make sure we have a mask for valid acceptor locations.
% Mask should be in a masks/ subdirectory and have the same filename as the
% original image stack with the addition of '_acceptorMask.tif'.
[path,filename,ext] = fileparts(fluorescenceImageStackFilePath);
acceptorMaskFilePath = fullfile(path,'masks',[filename '_acceptorMask.tif']);
if ~exist(acceptorMaskFilePath,'file')
    errordlg('ERROR: findAndProjectRoisForZmwAlexFret requires a valid TIFF image mask for acceptor locations.');
    return
end
acceptorMask = imread(acceptorMaskFilePath);

%% Show acceptor mask.
figure(roisFig);
subplot(2,4,6); cla; imshow(acceptorMask,[]); title('Acceptor Mask');

%% Bright field image registration.
% If masks for bright field image exist, use them to estimate ZMW locations
% and 2D affine transformation between cameras. Masks should be in a masks/
% subdirectory and have the same filename as the original bright field
% image with the addition of '_acceptorMask.tif' or '_donorMask.tif'.
if ~exist('brightFieldImageFilePath','var')
    [file,path] = uigetfile('*.tif','Load dual camera bright field image.');
    if ~isequal(file,0)
        brightFieldImageFilePath = [path file];
    end
end
if exist('brightFieldImageFilePath','var') && ~isempty(brightFieldImageFilePath)
    brightFieldImage = imread(brightFieldImageFilePath);
    % Metamorph stores dual camera images in one big concatenated image.
    % Acceptor 512x512 on left and donor 512x512 on right.
    donorBrightFieldImage = brightFieldImage(:,513:end);
    acceptorBrightFieldImage = brightFieldImage(:,1:512);
    % Show bright field images.
    figure(roisFig);
    subplot(2,4,1); cla; imshow(acceptorBrightFieldImage,[]); title('Acceptor BF');
    subplot(2,4,4); cla; imshow(donorBrightFieldImage,[]); title('Donor BF');
    % Check for masks.
    [path,filename,ext] = fileparts(brightFieldImageFilePath);
    donorBrightFieldMaskFilePath = fullfile(path,'masks',[filename '_donorMask.tif']);
    acceptorBrightFieldMaskFilePath = fullfile(path,'masks',[filename '_acceptorMask.tif']);
    if exist(donorBrightFieldMaskFilePath,'file') && exist(acceptorBrightFieldMaskFilePath,'file')
        donorBrightFieldMask = imread(donorBrightFieldMaskFilePath);
        acceptorBrightFieldMask = imread(acceptorBrightFieldMaskFilePath);
        % Find ZMW locations in masks.
        donorBrightFieldMaskRois = findRoisInMask(donorBrightFieldMask,minRoiArea,minRoiSeparation);
        acceptorBrightFieldMaskRois = findRoisInMask(acceptorBrightFieldMask,minRoiArea,minRoiSeparation);
        % Find acceptor to donor camera transform.
        donorZmwPts = vertcat(donorBrightFieldMaskRois.Centroid);
        acceptorZmwPts = vertcat(acceptorBrightFieldMaskRois.Centroid);
        opt.method = 'affine';
        opt.viz = 0;
        acceptorToDonorCameraTransform = cpd_register(donorZmwPts,acceptorZmwPts,opt);
        acceptorZmwPtsInDonorCamera = acceptorToDonorCameraTransform.Y;
    end
end

%% Show bright field images and masks.
if exist('acceptorBrightFieldImage','var') && ~isempty(acceptorBrightFieldImage)
    figure(roisFig);
    subplot(2,4,1); cla; imshow(acceptorBrightFieldImage,[]); title('Acceptor BF');
end
if exist('donorBrightFieldImage','var') && ~isempty(donorBrightFieldImage)
    figure(roisFig);
    subplot(2,4,4); cla; imshow(donorBrightFieldImage,[]); title('Donor BF');
end
if exist('acceptorBrightFieldMask','var') && ~isempty(acceptorBrightFieldMask)
    figure(roisFig);
    subplot(2,4,2); cla; imshow(acceptorBrightFieldMask,[]); title('Acceptor BF Mask');
end
if exist('donorBrightFieldMask','var') && ~isempty(donorBrightFieldMask)
    figure(roisFig);
    subplot(2,4,3); cla; imshow(donorBrightFieldMask,[]); title('Donor BF Mask');
end

%% Plot ROIs on top of bright field images and masks.
if exist('acceptorZmwPts','var')
    figure(roisFig);
    subplot(2,4,1); hold on; plot(acceptorZmwPts(:,1),acceptorZmwPts(:,2),'ro');
    subplot(2,4,2); hold on; plot(acceptorZmwPts(:,1),acceptorZmwPts(:,2),'ro');
end
if exist('donorZmwPts','var')
    figure(roisFig);
    subplot(2,4,3); hold on; plot(donorZmwPts(:,1),donorZmwPts(:,2),'gs');
    subplot(2,4,4); hold on; plot(donorZmwPts(:,1),donorZmwPts(:,2),'gs');
end
if exist('acceptorZmwPtsInDonorCamera','var')
    figure(roisFig);
    subplot(2,4,3); hold on; plot(acceptorZmwPtsInDonorCamera(:,1),acceptorZmwPtsInDonorCamera(:,2),'ro');
    subplot(2,4,4); hold on; plot(acceptorZmwPtsInDonorCamera(:,1),acceptorZmwPtsInDonorCamera(:,2),'ro');
end

%% Load image stack.
% Metamorph stores dual camera images in one big concatenated image.
% Acceptor 512x512 on left and donor 512x512 on right.
% Concatenated image stack will be 1024x512xN for N frames.
imageStack = loadTiffStack(fluorescenceImageStackFilePath,'Load fluorescene TIFF image stack.');
numFrames = int32(floor(size(imageStack,3)/2));

%% Stack averages for visualization.
acceptorInitialFrameAverage = double(squeeze(mean(imageStack(:,1:512,2:2:min([40,2*numFrames])),3)));
donorFrameAverage = double(squeeze(mean(imageStack(:,513:end,1:2:2*numFrames),3)));

%% Show acceptor and donor stack averages.
if exist('acceptorInitialFrameAverage','var') && ~isempty(acceptorInitialFrameAverage)
    figure(roisFig);
    subplot(2,4,5); cla; imshow(acceptorInitialFrameAverage,[]); title('Acceptor');
end
if exist('donorFrameAverage','var') && ~isempty(donorFrameAverage)
    figure(roisFig);
    subplot(2,4,8); cla; imshow(donorFrameAverage,[]); title('Donor');
end

%% Find valid acceptor molecules from mask.
acceptorRois = findRoisInMask(acceptorMask,minRoiArea,minRoiSeparation);
% Make sure acceptors are reasonably close to ZMW locations identified in
% the bright field image. If multiple acceptors are next to a ZMW, only
% keep the one that is closest.
if exist('acceptorZmwPts','var')
    acceptorPts = vertcat(acceptorRois.Centroid);
    acceptorNearestZmwIndices = zeros(size(acceptorPts,1),1);
    acceptorNearestZmwDistances = zeros(size(acceptorPts,1),1);
    for i = 1:size(acceptorPts,1)
        x0 = acceptorPts(i,1);
        y0 = acceptorPts(i,2);
        jmin = 0;
        d2min = 0;
        for j = 1:size(acceptorZmwPts,1)
            x1 = acceptorZmwPts(j,1);
            y1 = acceptorZmwPts(j,2);
            d2 = (x1-x0)^2+(y1-y0)^2;
            if (d2min == 0) || (d2 < d2min)
                jmin = j;
                d2min = d2;
            end
        end
        acceptorNearestZmwIndices(i) = jmin;
        acceptorNearestZmwDistances(i) = d2min;
    end
    for j = 1:size(acceptorZmwPts,1)
        nearbyAcceptors = find(acceptorNearestZmwIndices == j);
        if length(nearbyAcceptors) > 1
            [~,idx] = min(acceptorNearestZmwDistances(nearbyAcceptors));
            nearbyAcceptors(idx) = [];
            acceptorRois(nearbyAcceptors) = [];
            acceptorNearestZmwIndices(nearbyAcceptors) = [];
            acceptorNearestZmwDistances(nearbyAcceptors) = [];
        end
    end
end
% Refine acceptor locations with local Gaussian fits.
if acceptorRoiLocalSearch
    % !!! YOU MAY NEED TO ADJUST SIZE OF SEARCH BOX.
    roiSearchBox = [3,3];
    searchInEllipse = 0;
    visualize = 0;
    uiProgressBarMessage = 'Local refinement of acceptor ROI centers.';
    acceptorRois = findRoiGaussianCentersInImage(...
        acceptorInitialFrameAverage,acceptorRois,roiSearchBox,searchInEllipse,visualize,uiProgressBarMessage);
else
    for i = 1:length(acceptorRois)
        acceptorRois(i).Center = acceptorRois(i).Centroid;
    end
end

%% Plot ROIs on top of images and masks.
acceptorPts = vertcat(acceptorRois.Center);
figure(roisFig);
subplot(2,4,5); hold on;
    if exist('acceptorZmwPts','var')
        plot(acceptorZmwPts(:,1),acceptorZmwPts(:,2),'o','color',[.5,0,0]);
    end
    plot(acceptorPts(:,1),acceptorPts(:,2),'r+');
subplot(2,4,6); hold on;
    if exist('acceptorZmwPts','var')
        plot(acceptorZmwPts(:,1),acceptorZmwPts(:,2),'o','color',[.5,0,0]);
    end
    plot(acceptorPts(:,1),acceptorPts(:,2),'r+');

%% Find donor locations corresponding to acceptor locations.
if exist('acceptorToDonorCameraTransform','var')
    % Use bright field camera transform to guess initial donor locations.
    acceptorPts = vertcat(acceptorRois.Center);
    donorPts = cpd_transform(acceptorPts,acceptorToDonorCameraTransform);
    donorRois = acceptorRois;
    for i = 1:length(donorRois)
        donorRois(i).Center = donorPts(i,:);
    end
    % Remove ROIs whose donor locations are transformed outside the image
    % boundaries.
    donorPts = vertcat(donorRois.Center);
    outsideImageX = union(find(donorPts(:,1) < 0),find(donorPts(:,1) > 512));
    outsideImageY = union(find(donorPts(:,2) < 0),find(donorPts(:,2) > 512));
    outsideImage = union(outsideImageX,outsideImageY);
    if ~isempty(outsideImage)
        donorRois(outsideImage) = [];
        acceptorRois(outsideImage) = [];
    end
    % Refine donor locations with local Gaussian fits.
    if donorRoiLocalSearch
        % !!! YOU MAY NEED TO ADJUST SIZE OF SEARCH BOX.
        roiSearchBox = [3,3];
        searchInEllipse = 0;
        visualize = 0;
        uiProgressBarMessage = 'Local refinement of donor ROI centers.';
        donorRois = findRoiGaussianCentersInImage(...
            donorFrameAverage,donorRois,roiSearchBox,searchInEllipse,visualize,uiProgressBarMessage);
    end
else
    % Blindly search for donor locations near acceptor locations.
    % !!! YOU MAY NEED TO ADJUST SIZE OF SEARCH BOX.
    % !!! SEARCH BOX SHOULD BE BIG ENOUGH TO FIND TRANSFORMED DONOR
    % LOCATION, BUT NOT SO BIG SO AS TO INCLUDE NEIGHBORING ZMW LOCATIONS.
    roiSearchBox = [4,8];
    searchInEllipse = 0;
    visualize = 0;
    uiProgressBarMessage = 'Searching for donor ROIs.';
    donorRois = findRoiGaussianCentersInImage(...
        donorFrameAverage,acceptorRois,roiSearchBox,searchInEllipse,visualize,uiProgressBarMessage);
end
% Remove ROIs whose donor locations are transformed outside the image
% boundaries.
donorPts = vertcat(donorRois.Center);
outsideImageX = union(find(donorPts(:,1) < 0),find(donorPts(:,1) > 512));
outsideImageY = union(find(donorPts(:,2) < 0),find(donorPts(:,2) > 512));
outsideImage = union(outsideImageX,outsideImageY);
if ~isempty(outsideImage)
    donorRois(outsideImage) = [];
    acceptorRois(outsideImage) = [];
end

%% Plot ROIs on top of images and masks.
donorPts = vertcat(donorRois.Center);
acceptorPts = vertcat(acceptorRois.Center);
figure(roisFig);
subplot(2,4,5); cla; imshow(acceptorInitialFrameAverage,[]); title('Acceptor');
    hold on;
    if exist('acceptorZmwPts','var')
        plot(acceptorZmwPts(:,1),acceptorZmwPts(:,2),'o','color',[.5,0,0]);
    end
    plot(acceptorPts(:,1),acceptorPts(:,2),'r+');
subplot(2,4,6); cla; imshow(acceptorMask,[]); title('Acceptor Mask');
    hold on;
    if exist('acceptorZmwPts','var')
        plot(acceptorZmwPts(:,1),acceptorZmwPts(:,2),'o','color',[.5,0,0]);
    end
    plot(acceptorPts(:,1),acceptorPts(:,2),'r+');
subplot(2,4,8); cla; imshow(donorFrameAverage,[]); title('Donor');
    hold on;
    if exist('donorZmwPts','var')
        plot(donorZmwPts(:,1),donorZmwPts(:,2),'s','color',[0,.5,0]);
    end
    plot(donorPts(:,1),donorPts(:,2),'g+');
subplot(2,4,7); cla; imshow(acceptorInitialFrameAverage,[]); title('Acceptor (w/ Donor pts)');
    hold on; plot(donorPts(:,1),donorPts(:,2),'g+');
    plot(acceptorPts(:,1),acceptorPts(:,2),'r+');
    % Lines between donor-acceptor pairs.
    x = nan(3*size(donorPts,1),1);
    y = nan(3*size(donorPts,1),1);
    x(1:3:end) = donorPts(:,1);
    x(2:3:end) = acceptorPts(:,1);
    y(1:3:end) = donorPts(:,2);
    y(2:3:end) = acceptorPts(:,2);
    plot(x,y,'y-');

%% Project ROIs.
rois = repmat(struct(),size(acceptorRois));
wb = waitbar(0,'Projecting ROIs.');
numFrames = int32(floor(size(imageStack,3)/2));
for i = 1:length(rois)
    rois(i).acceptorCenter = acceptorRois(i).Center;
    rois(i).donorCenter = donorRois(i).Center;
    rois(i).acceptorPixels = findRoiPixels(rois(i).acceptorCenter,roiRadius,512,512);
    rois(i).donorPixels = findRoiPixels(rois(i).donorCenter,roiRadius,512,512);
    rois(i).zprojDD = zeros(numFrames,1);
    rois(i).zprojDA = zeros(numFrames,1);
    rois(i).zprojAA = zeros(numFrames,1);
    rois(i).zprojAD = zeros(numFrames,1);
    numAcceptorPixels = size(rois(i).acceptorPixels,1);
    if numAcceptorPixels
        for j = 1:numAcceptorPixels
            x = rois(i).acceptorPixels(j,1);
            y = rois(i).acceptorPixels(j,2);
            row = y;
            col = x;
            rois(i).zprojDA = rois(i).zprojDA + double(squeeze(imageStack(row,col,1:2:2*numFrames)));
            rois(i).zprojAA = rois(i).zprojAA + double(squeeze(imageStack(row,col,2:2:2*numFrames)));
        end
        rois(i).zprojDA = rois(i).zprojDA./numAcceptorPixels;
        rois(i).zprojAA = rois(i).zprojAA./numAcceptorPixels;
    end
    numDonorPixels = size(rois(i).donorPixels,1);
    if numDonorPixels
        for j = 1:numDonorPixels
            x = rois(i).donorPixels(j,1);
            y = rois(i).donorPixels(j,2);
            row = y;
            col = x+512;
            rois(i).zprojDD = rois(i).zprojDD + double(squeeze(imageStack(row,col,1:2:2*numFrames)));
            rois(i).zprojAD = rois(i).zprojAD + double(squeeze(imageStack(row,col,2:2:2*numFrames)));
        end
        rois(i).zprojDD = rois(i).zprojDD./numDonorPixels;
        rois(i).zprojAD = rois(i).zprojAD./numDonorPixels;
    end
    if mod(i,10) == 0
        waitbar(double(i)/length(rois),wb);
    end
end
close(wb);

%% Misc stuff.
info = parseFileName(fluorescenceImageStackFilePath,{'fcAMP'});
for i = 1:length(rois)
    rois(i).info = info;
    rois(i).isSelected = 0;
    rois(i).isFlagged = 0;
    rois(i).bleachFramesAA = [];
    rois(i).idealDA = [];
    rois(i).idealAA = [];
end

%% Save ROIs.
disp('Saving ROIs...');
[path,filename,ext] = fileparts(fluorescenceImageStackFilePath);
roisDir = fullfile(path,'rois');
if ~exist(roisDir,'dir')
    mkdir(roisDir);
end
save(fullfile(roisDir,[filename '.mat']),'rois');
disp('Finished saving ROIs.');

%% Save ROIs figure.
figsDir = fullfile(roisDir,'figs');
if ~exist(figsDir,'dir')
    mkdir(figsDir);
end
savefig(roisFig,fullfile(figsDir,[filename '.fig']));

end
