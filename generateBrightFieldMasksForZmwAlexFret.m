function generateBrightFieldMasksForZmwAlexFret(brightFieldImageFilePath)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com

if ~exist('brightFieldImageFilePath','var') || isempty(brightFieldImageFilePath)
    [file,path] = uigetfile('*.tif','Load bright field image.');
    if ~isequal(file,0)
        brightFieldImageFilePath = [path file];
    end
end

if exist('brightFieldImageFilePath','var') && ~isempty(brightFieldImageFilePath)
    brightFieldImage = imread(brightFieldImageFilePath);
    donorBrightFieldImage = brightFieldImage(:,513:end);
    acceptorBrightFieldImage = brightFieldImage(:,1:512);
    [path,filename,ext] = fileparts(brightFieldImageFilePath);
    maskDir = fullfile(path,'masks');
    if ~exist(maskDir,'dir')
        mkdir(maskDir);
    end
    imwrite(donorBrightFieldImage,fullfile(maskDir,[filename '_donorMask.tif']));
    imwrite(acceptorBrightFieldImage,fullfile(maskDir,[filename '_acceptorMask.tif']));
end

end
