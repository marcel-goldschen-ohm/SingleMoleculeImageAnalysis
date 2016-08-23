function info = parseFileName(pathToFile,ligands,verbose)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com
%
% Parse file path for specific parameters.
%
% args:
% -----
% pathToFile:   Path to file to parse.
% ligands:      Cell array of ligand names to search for in file name.
%
% return:
% -------
% info:         Map container with parsed fields.

%% Init.
info = containers.Map;
if ~exist('pathToFile','var') || isempty(pathToFile)
    % Get filename if it was not already given.
    [file,path] = uigetfile('*.tif','Select file.');
    if isequal(file,0)
        errordlg('ERROR: parseFileName requires a valid file.');
        return
    end
    pathToFile = [path file];
else
    seps = strfind(pathToFile,filesep);
    if isempty(seps)
        path = '';
        file = pathToFile;
    else
        path = pathToFile(1:seps(end));
        file = pathToFile(seps(end)+1:end);
    end
end
if ~exist('ligands','var')
    ligands = {};
end
if ~exist('verbose','var')
    verbose = 0;
end

%% Date (e.g. '2016-01-11')
% Could be in file path, but not file name.
info('Date') = '';
pos = regexp(pathToFile,'(\d\d\d\d)-(\d\d)-(\d\d)');
if ~isempty(pos)
    if length(pos) > 1
        pos = pos(end);
    end
    info('Date') = pathToFile(pos:pos+9);
end
if verbose
    disp(['Date: ' info('Date')]);
end

%% Image # (e.g. 'img003_')
info('Image #') = 0;
pos = regexp(file,'img(\d\d\d)_');
if ~isempty(pos)
    if length(pos) > 1
        pos = pos(1);
    end
    info('Image #') = uint32(str2num(file(pos+3:pos+5)));
end
if verbose
    disp(['Image #: ' num2str(info('Image #'))]);
end

%% Location # (e.g. '_loc002_')
info('Location #') = 0;
pos = regexp(file,'_loc(\d\d\d)_');
if ~isempty(pos)
    if length(pos) > 1
        pos = pos(1);
    end
    info('Location #') = uint32(str2num(file(pos+4:pos+6)));
end
if verbose
    disp(['Location #: ' num2str(info('Location #'))]);
end

%% ZMW diameter (e.g. '_150nm_' or '_d150nm_')
info('ZMW diameter (nm)') = 0;
pos = regexp(file,'_(d|ZMW)?(\d?\d\d)nm_');
if ~isempty(pos)
    if length(pos) > 1
        pos = pos(1);
    end
    if file(pos+1) == 'd'
        info('ZMW diameter (nm)') = uint32(str2num(file(pos+2:pos+4)));
    elseif strcmp(file(pos+1:pos+3),'ZMW')
        info('ZMW diameter (nm)') = uint32(str2num(file(pos+4:pos+6)));
    else
        info('ZMW diameter (nm)') = uint32(str2num(file(pos+1:pos+3)));
    end
end
if verbose
    disp(['ZMW diameter (nm): ' num2str(info('ZMW diameter (nm)'))]);
end

%% Excitation and frame rate (e.g. '_532nm45mW100ms_' or for ALEX: '_532nm45mW100ms640nm35mW100ms_')
info('Excitation Wavelength (nm)') = [];
info('Excitation Power (mW)') = [];
info('Excitation Exposure Time (ms)') = [];
pos = regexp(file,'([0-9]+)nm([0-9]+)mW([0-9]+)(ms|s)');
if ~isempty(pos)
    for k = 1:length(pos)
        pos1 = pos(k);
        pos2 = strfind(file(pos1:end),'nm');
        if length(pos2) > 1
            pos2 = pos2(1);
        end
        pos2 = pos1-1+pos2;
        pos3 = strfind(file(pos2:end),'mW');
        if length(pos3) > 1
            pos3 = pos3(1);
        end
        pos3 = pos2-1+pos3;
        pos4 = regexp(file(pos3:end),'(ms|s)_');
        if length(pos4) > 1
            pos4 = pos4(1);
        end
        pos4 = pos3-1+pos4;
        info('Excitation Wavelength (nm)') = ...
            [info('Excitation Wavelength (nm)'),str2num(file(pos1:pos2-1))];
        info('Excitation Power (mW)') = ...
            [info('Excitation Power (mW)'),str2num(file(pos2+2:pos3-1))];
        info('Excitation Exposure Time (ms)') = ...
            [info('Excitation Exposure Time (ms)'),str2num(file(pos3+2:pos4-1))];
        if strcmp(file(pos4:pos4+1),'ms')
            info('Excitation Exposure Time (ms)') = ...
                [info('Excitation Exposure Time (ms)'),str2num(file(pos3+2:pos4-1))];
        elseif file(pos4) == 's'
            info('Excitation Exposure Time (ms)') = ...
                [info('Excitation Exposure Time (ms)'),str2num(file(pos3+2:pos4-1))*1000];
        end
    end
end
if verbose
    disp(['Excitation Wavelength (nm): ' num2str(info('Excitation Wavelength (nm)'))]);
    disp(['Excitation Power (mW): ' num2str(info('Excitation Power (mW)'))]);
    disp(['Excitation Exposure Time (ms): ' num2str(info('Excitation Exposure Time (ms)'))]);
end

%% Ligand concentration (e.g. '_10uM_fcAMP')
for k = 1:length(ligands)
    ligand = ligands{k};
    info([ligand ' (M)']) = 0;
    pos = regexp(file,['_([0-9]+)[fpnum]?M_' ligand]);
    if ~isempty(pos)
        if length(pos) > 1
            pos = pos(1);
        end
        pos2 = pos-1+regexp(file(pos:end),['[fpnum]?M_' ligand]);
        info([ligand ' (M)']) = str2num(file(pos+1:pos2-1));
        if file(pos2) == 'f'
            info([ligand ' (M)']) = info([ligand ' (M)'])*1e-15;
        elseif file(pos2) == 'p'
            info([ligand ' (M)']) = info([ligand ' (M)'])*1e-12;
        elseif file(pos2) == 'n'
            info([ligand ' (M)']) = info([ligand ' (M)'])*1e-9;
        elseif file(pos2) == 'u'
            info([ligand ' (M)']) = info([ligand ' (M)'])*1e-6;
        elseif file(pos2) == 'm'
            info([ligand ' (M)']) = info([ligand ' (M)'])*1e-3;
        end
    end
    if verbose
        disp([ligand ' (M): ' num2str(info([ligand ' (M)']))]);
    end
end

end
