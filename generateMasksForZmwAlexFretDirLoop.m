   function generateMasksForZmwAlexFretDirLoop(dirpath)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com

if ~exist('dirpath','var') || isempty(dirpath)
    dirpath = uigetdir();
end
files = dir(dirpath);
for i = 1:length(files)
    if ~files(i).isdir
        [path,filename,ext] = fileparts(files(i).name);
        if ~strcmp(filename(1),'.') && strcmp(ext,'.tif') && isempty(strfind(filename,'-file'))
            if ~isempty(strfind(filename,'_BF'))
                generateBrightFieldMasksForZmwAlexFret(fullfile(dirpath,[filename ext]));
            else
                generateAcceptorMaskForZmwAlexFret(fullfile(dirpath,[filename ext]));
            end
        end
    end
end

end
