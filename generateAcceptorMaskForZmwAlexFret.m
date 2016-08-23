function generateAcceptorMaskForZmwAlexFret(fluorescenceImageStackFilePath)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com

if ~exist('fluorescenceImageStackFilePath','var') || isempty(fluorescenceImageStackFilePath)
    [file,path] = uigetfile('*.tif','Load image stack.');
    if ~isequal(file,0)
        fluorescenceImageStackFilePath = [path file];
    end
end

if exist('fluorescenceImageStackFilePath','var') && ~isempty(fluorescenceImageStackFilePath)
    imageStack = loadTiffStack(fluorescenceImageStackFilePath,'',40);
    acceptorInitialFrameAverage = uint16(squeeze(mean(imageStack(:,1:512,2:2:end),3)));
    [path,filename,ext] = fileparts(fluorescenceImageStackFilePath);
    maskDir = fullfile(path,'masks');
    if ~exist(maskDir,'dir')
        mkdir(maskDir);
    end
    imwrite(acceptorInitialFrameAverage,fullfile(maskDir,[filename '_acceptorMask.tif']));
end

end
