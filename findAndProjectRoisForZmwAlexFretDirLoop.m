function findAndProjectRoisForZmwAlexFretDirLoop(...
    fluorescenceImageStacksDirPath,brightFieldImagesDirPath,...
    minRoiArea,minRoiSeparation,acceptorRoiLocalSearch,donorRoiLocalSearch,roiRadius)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com

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

%% Loop through files in directories and find and project ROIs.
if ~exist('fluorescenceImageStacksDirPath','var') || isempty(fluorescenceImageStacksDirPath)
    fluorescenceImageStacksDirPath = uigetdir(pwd,'Select directory with fluorescence image stacks.');
end
if ~exist('brightFieldImagesDirPath','var') || isempty(brightFieldImagesDirPath)
    brightFieldImagesDirPath = uigetdir(pwd,'Select directory with bright field images.');
end
if isequal(brightFieldImagesDirPath,0)
    brightFieldImagesDirPath = fluorescenceImageStacksDirPath;
end
fluorescenceImageStackFiles = dir(fluorescenceImageStacksDirPath);
brightFieldImageFiles = dir(brightFieldImagesDirPath);
for i = 1:length(fluorescenceImageStackFiles)
    if ~fluorescenceImageStackFiles(i).isdir
        [path,filename,ext] = fileparts(fluorescenceImageStackFiles(i).name);
        if ~strcmp(filename(1),'.') && strcmp(ext,'.tif') && isempty(strfind(filename,'_BF')) && isempty(strfind(filename,'-file'))
            fluorescenceImageStackFilePath = fullfile(fluorescenceImageStacksDirPath,fluorescenceImageStackFiles(i).name);
            stackInfo = parseFileName(fluorescenceImageStackFilePath);
            for j = 1:length(brightFieldImageFiles)
                if ~brightFieldImageFiles(j).isdir
                    [path,filename,ext] = fileparts(brightFieldImageFiles(j).name);
                    if ~strcmp(filename(1),'.') && strcmp(ext,'.tif') && ~isempty(strfind(filename,'_BF')) && isempty(strfind(filename,'-file'))
                        brightFieldImageFilePath = fullfile(brightFieldImagesDirPath,brightFieldImageFiles(j).name);
                        brightFieldInfo = parseFileName(brightFieldImageFilePath);
                        if stackInfo('Location #') == brightFieldInfo('Location #')
                            disp('----------');
                            disp(brightFieldImageFiles(j).name);
                            disp(fluorescenceImageStackFiles(i).name);
                            findAndProjectRoisForZmwAlexFret(...
                                fluorescenceImageStackFilePath,brightFieldImageFilePath,...
                                minRoiArea,minRoiSeparation,acceptorRoiLocalSearch,donorRoiLocalSearch,roiRadius);
                            close all
                            break
                        end
                    end
                end
            end
        end
    end
end

end
