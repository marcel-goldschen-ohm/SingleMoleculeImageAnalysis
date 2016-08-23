function rois = findRoiGaussianCentersInImage(img,rois,roiSearchBox,searchInEllipse,visualize,uiProgressBarMessage)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com
%
% Refine ROI locations by local fitting with 2d Gaussian.
%
% Requires fmgaussfit (http://www.mathworks.com/matlabcentral/fileexchange/41938-fit-2d-gaussian-with-optimization-toolbox).
%
% args:
% -----
% img:                  2d image to be searched for ROI locations.
% rois:                 Struct array of starting ROIs (e.g. as returned by regionprops).
% roiSearchBox:         [dx,dy] for search bbox where x in [x0-dx,x0+dx]
%                       and y in [y0-dy,y0+dy].
% searchInEllipse:      If nonzero, limit search to Gaussians centered within
%                       the ellipse bounded by roiSearchBox.
% visualize:            If nonzero, show plot with fits and waitforbuttonpress
%                       between ROIs.
% uiProgressBarMessage: Message displayed in waitbar.
%
% return:
% -------
% rois:                 Updated struct array of ROIs containing results of
%                       Gaussian fits.

%% Init.
imageWidth = size(img,1);
imageHeight = size(img,2);
dx = roiSearchBox(1);
dy = roiSearchBox(2);
if ~exist('searchInEllipse','var')
    searchInEllipse = 0;
end
if ~exist('visualize','var')
    visualize = 0;
end
if ~exist('uiProgressBarMessage','var')
    uiProgressBarMessage = 'Finding ROI Gaussian centers.';
end

%% Make sure ROIs have either a Center, GaussianCenter or Centroid field.
hasCenter = any(strcmp('Center',fieldnames(rois)));
if ~hasCenter
    hasCentroid = any(strcmp('Centroid',fieldnames(rois)));
    if ~hasCentroid
        errordlg('ERROR: findRoisInImage requires ROIs with either a Center or Centroid field predefined.');
        return
    end
end

%% Find ROI 2D Gaussian centers.
if visualize
    h = figure;
    pos = get(h,'Position');
    pos(3) = 3*pos(4); % width:height == 3:1
    set(h,'Position',pos);
end
wb = waitbar(0,uiProgressBarMessage);
for i = 1:length(rois)
    if hasCenter
        x0 = rois(i).Center(1);
        y0 = rois(i).Center(2);
    else
        x0 = rois(i).Centroid(1);
        y0 = rois(i).Centroid(2);
    end
    % Set bounding box in which to search for ROI Gaussian.
    col0 = max([1,floor(x0-dx)]); % First column.
    row0 = max([1,floor(y0-dy)]); % First row.
    col1 = min([ceil(x0+dx),imageWidth]); % Last column.
    row1 = min([ceil(y0+dy),imageHeight]); % Last row.
    nrows = row1-row0+1; % Height in rows.
    ncols = col1-col0+1; % Width in columns.
    roiSubImage = img(row0:row1,col0:col1); % Sub-image to search for ROI.
    % Convert roisImage (x0,y0) to roiSubImage (roiX0,roiY0).
    roiX0 = x0-col0+1;
    roiY0 = y0-row0+1;
    % Plot ROI cmap with initial location marked.
    if visualize
        figure(h); subplot(1,3,1); cla; hold on;
        imagesc(roiSubImage);
        plot(roiX0,roiY0,'r+');
    end
    % Fit 2D gaussian to ROI.
    % fmgaussfit downloaded from:
    % http://www.mathworks.com/matlabcentral/fileexchange/41938-fit-2d-gaussian-with-optimization-toolbox
    [xx,yy] = meshgrid(1:size(roiSubImage,2),1:size(roiSubImage,1));
    warning('off','all');
    [fitresult,zfit,fiterr,zerr,resnorm,rr] = fmgaussfit(xx,yy,roiSubImage);
    warning('on','all');
    % fitresult = [A, ang, sx, sy, x0, y0, C]
    % gaussian2D = C + A * exp(-(( (x-x0).*cosd(ang)+(y-y0).*sind(ang))./sx).^2
    %                          -((-(x-x0).*sind(ang)+(y-y0).*cosd(ang))./sy).^2)
    A = fitresult(1);
    ang = fitresult(2);
    sx = fitresult(3);
    sy = fitresult(4);
    roiX1 = fitresult(5);
    roiY1 = fitresult(6);
    C = fitresult(7);
    % Convert roiSubImage (roiX1,roiY1) back to roisImage (x1,y1).
    x1 = roiX1+col0-1;
    y1 = roiY1+row0-1;
    success = 0;
    if searchInEllipse
        % Ellipse bounded by search box is ((x-x0)/a)^2 +((y-y0)/b)^2 = 1
        % If test point (x,y) in above equation < 1, it's inside the ellipse.
        % If test point (x,y) in above equation > 1, it's outside the ellipse.
        e = ((x1-x0)/(ncols/2))^2 +((y1-y0)/(nrows/2))^2;
        if e < 1
            success = 1;
        end
    else
        % Make sure center (x1,y1) is within the search box.
        if (x1 >= col0) && (x1 <= col1) && (y1 >= row0) && (y1 <= row1)
            success = 1;
        end
    end
    if success
        rois(i).Center = [x1,y1];
        rois(i).GaussianSigma = [sx,sy];
        rois(i).GaussianAmplitude = A;
        rois(i).GaussianBackground = C;
        rois(i).GaussianAngle = ang;
        if visualize
            figure(h); subplot(1,3,1); hold on;
            plot(roiX1,roiY1,'ro');
        end
    else
        % Keep original ROI location and zero other parameters.
        rois(i).Center = [x0,y0];
        rois(i).GaussianSigma = [0,0];
        rois(i).GaussianAmplitude = 0;
        rois(i).GaussianBackground = 0;
        rois(i).GaussianAngle = 0;
    end
    % Show fit over ROI heightmap.
    if visualize
        figure(h); subplot(1,3,2); cla; hold on;
        surf(xx,yy,roiSubImage);
        surf(xx,yy,zfit,'EdgeColor',[1,0,0],'FaceColor','none');
        % Symmetric gaussian.
        zz = C+A.*exp(-((xx-roiX1)./sigma).^2-((yy-roiY1)./sigma).^2);
        subplot(1,3,3); cla; hold on;
        surf(xx,yy,roiSubImage);
        surf(xx,yy,zz,'EdgeColor',[1,0,0],'FaceColor','none');
    end
    % Increment waitbar.
    waitbar(double(i)/length(rois),wb);
    if visualize
        waitforbuttonpress
    end
end
close(wb);

end
