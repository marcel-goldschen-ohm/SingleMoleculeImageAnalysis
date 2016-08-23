function [imageStack,info] = loadTiffStack(pathToFile,uiLoadMesssage,maxFrames)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com
%
% Load TIFF image stack from file.
%
% args:
% -----
% pathToFile:       Optional file path of TIFF image stack to be loaded.
%                   If not exist or empty --> popup a file selection UI.
% uiLoadMesssage:	Optional message displayed in file selection UI.
%                   If not exist or empty --> use default message.
%
% return:
% -------
% imageStack:       (width x height x frames) 3d array.
% info:             Struct array of info for each frame.

%% Init.
imageStack = [];
info = [];
defaultUiLoadMesssage = 'Load TIFF stack.';
if ~exist('maxFrames','var')
    maxFrames = 0;
end

%% Popup file selection dialog if file not given.
if ~exist('pathToFile','var') || isempty(pathToFile)
    % Get filename if it was not already given.
    if ~exist('loadMesssage','var') || isempty(uiLoadMesssage)
        uiLoadMesssage = defaultUiLoadMesssage;
    end
    [file,path] = uigetfile('*.tif',uiLoadMesssage);
    if isequal(file,0)
        errordlg('ERROR: loadTiffStack requires a valid TIFF file.');
        return
    end
    pathToFile = [path file];
    fileName = file;
else
    seps = strfind(pathToFile,filesep);
    fileName = pathToFile(seps(end)+1:end);
end

%% File label (displayed in waitbar).
% Convert '_' to ' '.
fileLabel = fileName;
idx = strfind(fileName,'_');
for k = 1:length(idx)
    fileLabel(idx(k)) = ' ';
end

%% Start timer.
tic

%% Get some info about TIFF file.
info = imfinfo(pathToFile);
numFrames = numel(info);
if numFrames == 0
    errordlg('ERROR: Zero frames detected in image stack.');
    return
end
if maxFrames > 0 && maxFrames < numFrames
    numFrames = maxFrames;
end
width = info(1).Width;
height = info(1).Height;
bpp = info(1).BitDepth;
color = info(1).ColorType;
if isempty(find([8 16 32] == bpp,1))
    errordlg(['ERROR: Unsupported image bit depth: ' num2str(bpp) ' (only 8, 16 and 32 bpp supported)']);
    return
end
if ~strcmp(color,'grayscale')
    errordlg(['ERROR: Unsupported image format: ' color ' (only grayscale supported)']);
    return
end
fmt = ['uint' num2str(bpp)];
infoStr = [num2str(width) 'x' num2str(height) 'x' num2str(numFrames) '@' num2str(bpp)];
disp('Loading TIFF stack...');
disp(['--> File: ' pathToFile]);
disp(['--> Info: ' infoStr]);

%% Allocate memory for entire image stack.
rows = height;
cols = width;
imageStack = zeros(rows,cols,numFrames,fmt);

%% Load stack one frame at a time.
% !!! Updating waitbar is expensive, so do it sparingly.
wb = waitbar(0,[fileLabel ' (' infoStr ')']);
oneTenthOfTotalFrames = floor(double(numFrames)/10);
for frame = 1:numFrames
    imageStack(:,:,frame) = imread(pathToFile,frame,'Info',info);
    if mod(frame, oneTenthOfTotalFrames) == 0
        waitbar(double(frame)/numFrames,wb);
    end
end
close(wb);

%% Check for additional files belonging to this stack.
% e.g. myFile.tif, myFile-file002.tif, myFile-file003.tif, ...
done = 0;
if maxFrames > 0
    maxFrames = maxFrames - numFrames;
    if maxFrames <= 0
        done = 1;
    end
end
if ~done
    if strcmp(pathToFile(end-11:end-7),'-file')
        fileNum = str2num(pathToFile(end-6:end-4));
        pathToNextFile = [pathToFile(1:end-7) num2str(fileNum+1,'%03u') '.tif'];
    else
        pathToNextFile = [pathToFile(1:end-4) '-file002.tif'];
    end
    if exist(pathToNextFile,'file') == 2
        nextStack = loadTiffStack(pathToNextFile,'',maxFrames);
        imageStack = cat(3,imageStack,nextStack);
    end
end

%% Elapsed time.
toc
disp('... Finished loading TIFF stack.');

end
