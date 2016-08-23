function roiProjectionViewerForZmwAlexFret()
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com

%% Main UI layout.
mainFig = figure('units','normalized','position',[0.05, 0.05, 0.80, 0.70]);
set(mainFig,'units','pixels');
pos = get(mainFig,'position');
w = pos(3);
h = pos(4);
% | buttons/controls | 4x1 image stacks | 4x1 ROI projections |
btnw = 200; % Buttons/controls panel width (px).
axw = 40; % Plot axis width (px).
s = 15; % Spacer (px).
rwh = 6; % Projection plot width/height ratio.
w = floor((w-s-btnw-axw-axw-s)/(1+rwh)); % Available unit width (px).
h = floor((h-axw-3*s-axw)/4); % Available unit height (px).
d = min([w,h]); % Unit square dimension (px).
w = s+btnw+axw+d+axw+rwh*d+s; % Contents width for given unit square dimension (px).
h = axw+4*d+3*s+axw; % Contents height for given unit square dimension (px).
imgAxesDD = axes('units','pixels','position',[s+btnw+axw,axw+3*(d+s),d,d]);
imgAxesDA = axes('units','pixels','position',[s+btnw+axw,axw+2*(d+s),d,d]);
imgAxesAA = axes('units','pixels','position',[s+btnw+axw,axw+1*(d+s),d,d]);
imgAxesAD = axes('units','pixels','position',[s+btnw+axw,axw+0*(d+s),d,d]);
projAxesDD = axes('units','pixels','position',[s+btnw+axw+d+axw,axw+3*(d+s),rwh*d,d]);
projAxesDA = axes('units','pixels','position',[s+btnw+axw+d+axw,axw+2*(d+s),rwh*d,d]);
projAxesAA = axes('units','pixels','position',[s+btnw+axw+d+axw,axw+1*(d+s),rwh*d,d]);
projAxesAD = axes('units','pixels','position',[s+btnw+axw+d+axw,axw+0*(d+s),rwh*d,d]);
set(mainFig,'position',[pos(1),pos(2),w,h]);
linkaxes([projAxesDD,projAxesDA,projAxesAA,projAxesAD],'x');

%% Buttons & Controls.
lineh = 18;
y = h;
% ROI analysis i/o.
y = y-2*lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','load','horizontalAlignment','center','callback',{@loadRois},'position',[s,y,floor(btnw/3),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','append','horizontalAlignment','center','callback',{@appendRois},'position',[s+floor(btnw/3),y,floor(btnw/3),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','save','horizontalAlignment','center','callback',{@saveRois},'position',[s+2*floor(btnw/3),y,floor(btnw/3),lineh]);
% ROI traversal.
y = y-2*lineh;
roiIndexEdit = uicontrol('parent',mainFig,'style','edit','string','0','horizontalAlignment','left','callback',{@roiIndexChanged},'position',[s,y,floor(btnw/3),lineh]);
numRoisText = uicontrol('parent',mainFig,'style','text','string','/ 0 ROIs','horizontalAlignment','left','position',[s+floor(btnw/3),y,btnw-4*lineh-floor(btnw/3),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','<','callback',{@prevRoi},'position',[s+btnw-4*lineh,y,2*lineh,lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','>','callback',{@nextRoi},'position',[s+btnw-2*lineh,y,2*lineh,lineh]);
% ROI selection.
y = y-2*lineh;
isRoiSelectedCheckBox = uicontrol('parent',mainFig,'style','checkbox','string','ROI selected (0 total)','value',0,'horizontalAlignment','left','callback',{@roiOptionsChanged},'position',[s,y,btnw,lineh]);
y = y-lineh;
isRoiFlaggedCheckBox = uicontrol('parent',mainFig,'style','checkbox','string','ROI flagged (0 total)','value',0,'horizontalAlignment','left','callback',{@roiOptionsChanged},'position',[s,y,btnw,lineh]);
y = y-lineh;
showOnlySelectedRoisCheckBox = uicontrol('parent',mainFig,'style','checkbox','string','show ONLY selected ROIs','value',0,'horizontalAlignment','left','position',[s,y,btnw,lineh]);
y = y-lineh;
showOnlyFlaggedRoisCheckBox = uicontrol('parent',mainFig,'style','checkbox','string','show ONLY flagged ROIs','value',0,'horizontalAlignment','left','position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','select all','horizontalAlignment','center','callback',{@selectAllRois},'position',[s,y,floor(btnw/3),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','unselect all','horizontalAlignment','center','callback',{@unselectAllRois},'position',[s+floor(btnw/3),y,floor(btnw/3),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','toggle selected','horizontalAlignment','center','callback',{@toggleSelectedRois},'position',[s+2*floor(btnw/3),y,floor(btnw/3),lineh]);
% Single AA bleach step.
y = y-2*lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','<','horizontalAlignment','center','callback',{@prevSingleBleachStepForCurrentRoiAA},'position',[s,y,2*lineh,lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','Find single AA bleach step','horizontalAlignment','center','callback',{@findSingleBleachStepForCurrentRoiAA},'position',[s+2*lineh,y,btnw-4*lineh,lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','>','horizontalAlignment','center','callback',{@nextSingleBleachStepForCurrentRoiAA},'position',[s+btnw-2*lineh,y,2*lineh,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','Set single AA bleach step','horizontalAlignment','center','callback',{@setSingleBleachStepForCurrentRoiAA},'position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','Find single AA bleach step for all ROIs','horizontalAlignment','center','callback',{@findSingleBleachStepForAllRoisAA},'position',[s,y,btnw,lineh]);
% DD cross talk in DA.
y = y-2*lineh;
uicontrol('parent',mainFig,'style','text','string','DA X-Talk','horizontalAlignment','left','position',[s,y,floor(0.75*btnw),lineh]);
crossTalkFractionEdit = uicontrol('parent',mainFig,'style','edit','string','0.107','horizontalAlignment','left','callback',{@roiOptionsChanged},'position',[s+floor(0.75*btnw),y,floor(0.25*btnw),lineh]);
% DD spline fit.
y = y-2*lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','fit DD spline','horizontalAlignment','center','callback',{@fitSplineDD},'position',[s,y,floor(btnw/2),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','clear DD spline','horizontalAlignment','center','callback',{@clearSplineDD},'position',[s+floor(btnw/2),y,floor(btnw/2),lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','text','string','# DD spline ctrl pts','horizontalAlignment','left','position',[s,y,floor(0.75*btnw),lineh]);
numsplineCtrlPtsEditDD = uicontrol('parent',mainFig,'style','edit','string','0','horizontalAlignment','left','position',[s+floor(0.75*btnw),y,floor(0.25*btnw),lineh]);
% DA spline fit.
y = y-2*lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','fit DA spline','horizontalAlignment','center','callback',{@fitSplineDA},'position',[s,y,floor(btnw/2),lineh]);
uicontrol('parent',mainFig,'style','pushbutton','string','clear DA spline','horizontalAlignment','center','callback',{@clearSplineDA},'position',[s+floor(btnw/2),y,floor(btnw/2),lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','text','string','# DA spline ctrl pts','horizontalAlignment','left','position',[s,y,floor(0.75*btnw),lineh]);
numsplineCtrlPtsEditDA = uicontrol('parent',mainFig,'style','edit','string','0','horizontalAlignment','left','position',[s+floor(0.75*btnw),y,floor(0.25*btnw),lineh]);
% vbFRET i/o.
y = y-2*lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','export DA for vbFRET','horizontalAlignment','center','callback',{@exportCurrentRoiForVbFretDA},'position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','import DA from vbFRET','horizontalAlignment','center','callback',{@importIdealForCurrentRoiFromVbFretDA},'position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','export all selected DA for vbFRET','horizontalAlignment','center','callback',{@exportAllSelectedRoisForVbFretDA},'position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','import all selected DA from vbFRET','horizontalAlignment','center','callback',{@importIdealForAllSelectedRoisFromVbFretDA},'position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','export all flagged DA for vbFRET','horizontalAlignment','center','callback',{@exportAllFlaggedRoisForVbFretDA},'position',[s,y,btnw,lineh]);
y = y-lineh;
uicontrol('parent',mainFig,'style','pushbutton','string','import all selected DA from vbFRET','horizontalAlignment','center','callback',{@importIdealForAllFlaggedRoisFromVbFretDA},'position',[s,y,btnw,lineh]);

%% Handle keyboard input.
% !!! This does NOT work while zoomed on any of the plots,
%     or indeed if any other uimode is active.
%     For a possible work around, see:
%       http://undocumentedmatlab.com/blog/enabling-user-callbacks-during-zoom-pan
set(mainFig,'KeyPressFcn',@keyPress);

    function keyPress(obj,event)
        key = event.Key; %disp(key);
        if strcmp(key,'leftarrow')
            prevRoi();
        elseif strcmp(key,'rightarrow')
            nextRoi();
        elseif strcmp(key,'uparrow')
            set(isRoiSelectedCheckBox,'value',1);
            roiOptionsChanged();
        elseif strcmp(key,'downarrow')
            set(isRoiSelectedCheckBox,'value',0);
            roiOptionsChanged();
        elseif strcmp(key,'d')
            fitSplineDA;
        elseif strcmp(key,'e')
            clearSplineDA;   
        elseif strcmp(key,'f')
            flagged = get(isRoiFlaggedCheckBox,'value');
            if flagged == 0
                flagged = 1;
            else
                flagged = 0;
            end
            set(isRoiFlaggedCheckBox,'value',flagged);
            roiOptionsChanged();
        end
    end

%% Init.
roisFilePath = '';
rois = []; % Struct array (e.g. similar to that returned by regionprops).
% rois(k).acceptorCenter = [x,y]
% rois(k).donorCenter = [x,y]
% rois(k).acceptorPixels = [[x,y];[x,y];...]
% rois(k).donorPixels = [[x,y];[x,y];...]
% rois(k).zprojDD = []
% rois(k).zprojDA = []
% rois(k).zprojAA = []
% rois(k).zprojAD = []
% rois(k).isSelected = 0 or 1
% rois(k).isFlagged = 0 or 1
% rois(k).bleachFramesAA = []
% rois(k).idealDA = []
% rois(k).idealAA = []

%% ROIs analysis i/o.
    function loadRois(obj,event,filepath)
        rois = [];
        if exist('filepath','var')
            appendRois(obj,event,filepath);
        else
            appendRois(obj,event);
        end
    end

    function appendRois(obj,event,filepath)
        if ~exist('filepath','var') || isempty(filepath)
            [file,path] = uigetfile('*.mat','Load ROIs analysis (*.mat).');
            if isequal(file,0)
                return
            end
            filepath = [path file];
        end
        disp('Loading ROIs analysis...');
        temp = load(filepath,'rois');
        newRois = checkDataSanity(temp.rois);
        firstNewRoiIndex = length(rois)+1;
        rois = [rois;newRois];
        roisFilePath = filepath;
        disp('... Finished loading ROIs analysis.');
        goToRoi(firstNewRoiIndex);
    end

    function saveRois(obj,event,filepath)
        if ~exist('filepath','var') || isempty(filepath)
            defaultSaveFilePath = roisFilePath;
            if isempty(defaultSaveFilePath)
                defaultSaveFilePath = '*.mat';
            end
            [file,path] = uiputfile(defaultSaveFilePath,'Save ROIs analysis (*.mat).');
            if isequal(file,0)
                return
            end
            filepath = [path file];
        end
        disp('Saving ROIs analysis...');
        save(filepath,'rois');
        roisFilePath = filepath;
        disp('... Finished saving ROIs analysis.');
    end

    function saneRois = checkDataSanity(insaneRois)
        if ~isfield(insaneRois,'zprojDD')
            for k = 1:length(insaneRois)
                insaneRois(k).zprojDD = [];
            end
        end
        if ~isfield(insaneRois,'zprojDA')
            for k = 1:length(insaneRois)
                insaneRois(k).zprojDA = [];
            end
        end
        if ~isfield(insaneRois,'zprojAD')
            for k = 1:length(insaneRois)
                insaneRois(k).zprojAD = [];
            end
        end
        if ~isfield(insaneRois,'zprojAA')
            for k = 1:length(insaneRois)
                insaneRois(k).zprojAA = [];
            end
        end
        if ~isfield(insaneRois,'isSelected')
            for k = 1:length(insaneRois)
                insaneRois(k).isSelected = 0;
            end
        end
        if ~isfield(insaneRois,'isFlagged')
            for k = 1:length(insaneRois)
                insaneRois(k).isFlagged = 0;
            end
        end
        if ~isfield(insaneRois,'bleachFramesAA')
            for k = 1:length(insaneRois)
                insaneRois(k).bleachFramesAA = [];
            end
        end
        if ~isfield(insaneRois,'idealDA')
            for k = 1:length(insaneRois)
                insaneRois(k).idealDA = [];
            end
        end
        if ~isfield(insaneRois,'idealAA')
            for k = 1:length(insaneRois)
                insaneRois(k).idealAA = [];
            end
        end
        if ~isfield(insaneRois,'splineDD')
            for k = 1:length(insaneRois)
                insaneRois(k).splineDD = [];
            end
        end
        if ~isfield(insaneRois,'splineDA')
            for k = 1:length(insaneRois)
                insaneRois(k).splineDA = [];
            end
        end
        saneRois = insaneRois;
    end

%% Go to specified ROI.
    function goToRoi(roiIndex)
        % If ROI index is not given, grab it from the UI.
        if ~exist('roiIndex','var')
            roiIndex = currentRoiIndex();
        end
        % Check validity of ROI index.
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return;
        end
        % Update UI.
        numSelectedRois = length(find(vertcat(rois.isSelected) > 0));
        numFlaggedRois = length(find(vertcat(rois.isFlagged) > 0));
        set(roiIndexEdit,'string',num2str(roiIndex));
        set(numRoisText,'string',['/ ' num2str(numRois) ' ROIs']);
        set(isRoiSelectedCheckBox,'string',['ROI selected (' num2str(numSelectedRois) ' total)']);
        set(isRoiSelectedCheckBox,'value',rois(roiIndex).isSelected);
        set(isRoiFlaggedCheckBox,'string',['ROI flagged (' num2str(numFlaggedRois) ' total)']);
        set(isRoiFlaggedCheckBox,'value',rois(roiIndex).isFlagged);
        % Draw ROI image stacks.
        % ...
        % Draw ROI projections.
        drawRoiProjections(roiIndex);
        % Return UI control to mainFig (for keyboard interaction).
        figure(mainFig);
    end

%% Get current ROI index from UI.
    function roiIndex = currentRoiIndex()
        roiIndex = floor(str2num(get(roiIndexEdit,'string')));
    end

%% Update ROI index.
    function roiIndexChanged(obj,event)
        roiIndex = currentRoiIndex();
        goToRoi(roiIndex);
    end

%% Go to previous ROI.
    function prevRoi(obj,event)
        roiIndex = currentRoiIndex()-1;
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return;
        end
        if get(showOnlySelectedRoisCheckBox,'value') && ~rois(roiIndex).isSelected
            prevSelected = find(vertcat(rois(1:roiIndex).isSelected) > 0);
            if isempty(prevSelected)
                return
            end
            roiIndex = prevSelected(end);
        end
        if get(showOnlyFlaggedRoisCheckBox,'value') && ~rois(roiIndex).isFlagged
            prevFlagged = find(vertcat(rois(1:roiIndex).isFlagged) > 0);
            if isempty(prevFlagged)
                return
            end
            roiIndex = prevFlagged(end);
        end
        goToRoi(roiIndex);
    end

%% Go to next ROI.
    function nextRoi(obj,event)
        prevRoiIndex = currentRoiIndex();
        roiIndex = prevRoiIndex+1;
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return;
        end
        if get(showOnlySelectedRoisCheckBox,'value') && ~rois(roiIndex).isSelected
            nextSelected = find(vertcat(rois(roiIndex:end).isSelected) > 0);
            if isempty(nextSelected)
                return
            end
            roiIndex = roiIndex-1+nextSelected(1);
        end
        if get(showOnlyFlaggedRoisCheckBox,'value') && ~rois(roiIndex).isFlagged
            nextFlagged = find(vertcat(rois(roiIndex:end).isFlagged) > 0);
            if isempty(nextFlagged)
                return
            end
            roiIndex = roiIndex-1+nextFlagged(1);
        end
        goToRoi(roiIndex);
    end

%% Update ROI options.
    function roiOptionsChanged(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        rois(roiIndex).isSelected = get(isRoiSelectedCheckBox,'value');
        rois(roiIndex).isFlagged = get(isRoiFlaggedCheckBox,'value');
        goToRoi(roiIndex);
    end

    function selectAllRois(obj,event)
        for k = 1:length(rois)
            rois(k).isSelected = 1;
        end
        goToRoi(); % Redraw.
    end

    function unselectAllRois(obj,event)
        for k = 1:length(rois)
            rois(k).isSelected = 0;
        end
        goToRoi(); % Redraw.
    end

    function toggleSelectedRois(obj,event)
        for k = 1:length(rois)
            if rois(k).isSelected == 0
                rois(k).isSelected = 1;
            else
                rois(k).isSelected = 0;
            end
        end
        goToRoi(); % Redraw.
    end

%% Draw ROI projections.
    function drawRoiProjections(roiIndex)
        roi = rois(roiIndex);
        % Get roi projections.
        yDD = roi.zprojDD;
        yDA = roi.zprojDA;
        yAA = roi.zprojAA;
        yAD = roi.zprojAD;
        iDA = roi.idealDA;
        iAA = roi.idealAA;
        sDD = roi.splineDD;
        sDA = roi.splineDA;
        if length(iAA) ~= length(yAA) && ~isempty(roi.bleachFramesAA)
            iAA = buildIdealFromBleachSteps(yAA,roi.bleachFramesAA);
        end
        frames = 1:length(yDA);
        if size(frames) ~= size(yDA)
            frames = frames';
        end
        % Subtract donor cross talk from FRET channel.
        xtalk = str2num(get(crossTalkFractionEdit,'string'));
        if xtalk > 0
            if ~isempty(sDD) && length(sDD) == length(yDD)
                yDA = yDA-((yDD-sDD+mean(yDD)).*xtalk);
            else
                yDA = yDA-(yDD.*xtalk);
            end
        end
        % Subtract DA spline.
        if ~isempty(sDA) && length(sDA) == length(yDA)
            yDA = yDA-sDA;
        end
        % Plot roi projections.
        axes(projAxesDD); cla;
        plot(frames,yDD,'-','color',[0,.75,0]);
        ylabel('DD (au)');
        axes(projAxesDA); cla;
        plot(frames,yDA,'-','color',[.95,0,0]);
        ylabel('DA (au)');
        axes(projAxesAA); cla;
        plot(frames,yAA,'-','color',[.85,0,.85]);
        ylabel('AA (au)');
        axes(projAxesAD); cla;
        plot(frames,yAD,'-','Color',[0,.75,.75]);
        ylabel('AD (au)');
        xlabel('Frame');
        % DD spline.
        if ~isempty(sDD) && length(sDD) == length(yDD)
            axes(projAxesDD); hold on;
            plot(frames,sDD,'-','color',[0,1,0]);
        end
        % DA spline.
%         if ~isempty(sDA) && length(sDA) == length(yDA)
%             axes(projAxesDA); hold on;
%             plot(frames,sDA,'-','color',[1,1,0]);
%         end
        % AA ideal bleach step(s).
        if ~isempty(iAA) && length(iAA) == length(yAA)
            axes(projAxesAA); hold on;
            plot(frames,iAA,'-','color',[0,0,1]);
        end
        % DA baseline.
        if ~isempty(roi.bleachFramesAA)
            baselineDA = mean(yDA(roi.bleachFramesAA(1):end));
            axes(projAxesDA); hold on;
            plot([frames(1),frames(end)],[baselineDA,baselineDA],'--','color',[0,0,1]);
        end
        % DA idealization. To not view ideal, == length(yDA) && 0
        if ~isempty(iDA) && length(iDA) == length(yDA)
            axes(projAxesDA); hold on;
            plot(frames,iDA,'-','color',[0,0,1]);
        end
        % Overlay DA artificial FRET on DA projection.
        % Yeah, duh. It's exactly the same.
        if 0
            axes(projAxesDA); hold on;
            yA = yDA-baselineDA;
            yD = -yA+max(yA);
            fret = scaleIdealBinaryToData(yA./(yA+yD),yDA);
            plot(frames,fret,'-','color',[0,0,1]);
        end
%         % Info.
%         axes(projAxesDD);
%         title(['fcAMP (uM) = ' num2str(roi.fcAMP_uM) ', ZMW (nm) = ' num2str(roi.ZMW_diameter_nm)]);
    end

%% Build ideal from bleach steps.
    function ideal = buildIdealFromBleachSteps(data,bleachFrames)
        ideal = zeros(size(data));
        frameStarts = [1,bleachFrames];
        for k=1:length(frameStarts)
            firstFrame = frameStarts(k);
            if k < length(frameStarts)
                lastFrame = frameStarts(k+1)-1;
            else
                lastFrame = length(data);
            end
            ideal(firstFrame:lastFrame) = mean(data(firstFrame:lastFrame));
        end
    end

%% Find single bleach step in AA.
    function findSingleBleachStepAA(roiIndex,frames)
        % If ROI index is not given, grab it from the UI.
        if ~exist('roiIndex','var')
            roiIndex = currentRoiIndex();
        end
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        yAA = rois(roiIndex).zprojAA;
        % If no frames were specified, use all frames.
        if ~exist('frames','var') || isempty(frames)
            frames = 1:length(yAA);
        else
            yAA = yAA(frames);
        end
        % Find minimal SSE for bleach step func at each frame with data.
        fit = zeros(size(yAA));
        cost = zeros(size(yAA));
        for i = 2:length(yAA)
            fit(1:i-1) = mean(yAA(1:i-1));
            fit(i:end) = mean(yAA(i:end));
            cost(i) = sum((fit-yAA).^2);
        end
        [~,bleachFrame] = min(cost(2:end));
        bleachFrameClipRange = 1+bleachFrame;
        bleachFrameFullRange = frames(bleachFrameClipRange);
        rois(roiIndex).bleachFramesAA = [bleachFrameFullRange];
        rois(roiIndex).idealAA = zeros(size(rois(roiIndex).zprojAA));
        rois(roiIndex).idealAA(1:bleachFrameFullRange-1) = mean(yAA(1:bleachFrameClipRange-1));
        rois(roiIndex).idealAA(bleachFrameFullRange:end) = mean(yAA(bleachFrameClipRange:end));
        goToRoi(roiIndex);
    end

    function findSingleBleachStepForCurrentRoiAA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        findSingleBleachStepAA(roiIndex);
    end

    function prevSingleBleachStepForCurrentRoiAA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        roi = rois(roiIndex);
        numFrames = length(roi.zprojAA);
        if ~isempty(roi.bleachFramesAA)
            bleachFrame = roi.bleachFramesAA(1);
            if (bleachFrame > 2) && (bleachFrame <= numFrames)
                findSingleBleachStepAA(roiIndex,1:bleachFrame-1);
            end
        end
    end

    function nextSingleBleachStepForCurrentRoiAA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        roi = rois(roiIndex);
        numFrames = length(roi.zprojAA);
        if ~isempty(roi.bleachFramesAA)
            bleachFrame = roi.bleachFramesAA(1);
            if (bleachFrame > 1) && (bleachFrame < numFrames)
                findSingleBleachStepAA(roiIndex,bleachFrame+1:numFrames);
            end
        end
    end

    function setSingleBleachStepForCurrentRoiAA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        numFrames = length(rois(roiIndex).zprojAA);
        prompt = {['Enter acceptor bleach frame (2-' num2str(numFrames) '):']};
        title = 'Acceptor Bleach Frame';
        numLines = 1;
        defaultAnswer = {num2str(numFrames)};
        answer = inputdlg(prompt,title,numLines,defaultAnswer);
        if ~isempty(answer)
            bleachFrame = str2num(answer{1});
            if (bleachFrame >= 2) && (bleachFrame <= numFrames)
                rois(roiIndex).bleachFramesAA = [bleachFrame];
                goToRoi(); % Redraw.
            end
        end
    end

    function findSingleBleachStepForAllRoisAA(obj,event)
        for roiIndex=1:length(rois)
            findSingleBleachStepAA(roiIndex);
        end
        goToRoi(1);
    end

%% DD spline fit.
    function fitSplineDD(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        yDD = rois(roiIndex).zprojDD;
        numFrames = length(yDD);
        frames = [1:numFrames]';
        numCtrlPts = str2num(get(numsplineCtrlPtsEditDD,'string'));
        if numCtrlPts == 0
            numCtrlPts = max([2,int32(floor(double(numFrames)/100))]);
        end
        pp = splinefit(frames,yDD,linspace(0,numFrames,int32(numCtrlPts)));
        rois(roiIndex).splineDD = ppval(pp,frames);
        goToRoi(roiIndex);
    end

    function clearSplineDD(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        rois(roiIndex).splineDD = [];
        goToRoi(roiIndex);
    end

%% DA spline fit.
    function fitSplineDA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        roi = rois(roiIndex);
        if isempty(roi.bleachFramesAA)
            return
        end
        yDA = roi.zprojDA;
        yDD = roi.zprojDD;
        sDD = roi.splineDD;
        % Subtract donor cross talk from FRET channel.
        xtalk = str2num(get(crossTalkFractionEdit,'string'));
        if xtalk > 0
            if ~isempty(sDD) && length(sDD) == length(yDD)
                yDA = yDA-((yDD-sDD+mean(yDD)).*xtalk);
            else
                yDA = yDA-(yDD.*xtalk);
            end
        end
        % Subtract non-constant baseline from FRET channel.
        amp = roi.baselineDecayAmpDA;
        if amp ~= 0
            tau = roi.baselineDecayTauDA;
            yDA = subtractBaselineDecay(yDA,amp,tau);
        end
        bleachFrame = roi.bleachFramesAA(1);
        numFrames = length(yDA);
        x = [bleachFrame+1:numFrames]';
        y = yDA(x);
        numCtrlPts = str2num(get(numsplineCtrlPtsEditDA,'string'));
        if numCtrlPts == 0
            numCtrlPts = max([2,int32(floor(double(length(y))/100))]);
        end
        pp = splinefit(x,y,linspace(x(1),x(end),int32(numCtrlPts)));
        ys = ppval(pp,x);
        m = ys(1)-ys(2);
        yR0 = max(yDA(1:bleachFrame))-min(yDA(1:bleachFrame));
        if abs(m*bleachFrame) < yR0
            rois(roiIndex).splineDA = zeros(size(yDA));
            rois(roiIndex).splineDA(x) = ys;
            for k = 1:bleachFrame
                rois(roiIndex).splineDA(k) = ys(1)+m*(bleachFrame-k+1);
            end
            goToRoi(roiIndex);
        end
    end

    function clearSplineDA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        rois(roiIndex).splineDA = [];
        goToRoi(roiIndex);
    end

%% vbFRET i/o.
    function exportCurrentRoiForVbFretDA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        exportForVbFretDA([roiIndex]);
    end

    function importIdealForCurrentRoiFromVbFretDA(obj,event)
        roiIndex = currentRoiIndex();
        numRois = length(rois);
        if (roiIndex <= 0) || (roiIndex > numRois)
            return
        end
        importFromVbFretDA([roiIndex]);
    end

    function exportAllSelectedRoisForVbFretDA(obj,event)
        roiIndices = find(vertcat(rois.isSelected) > 0);
        if ~isempty(roiIndices)
            exportForVbFretDA(roiIndices);
        end
    end

    function importIdealForAllSelectedRoisFromVbFretDA(obj,event)
        roiIndices = find(vertcat(rois.isSelected) > 0);
        if ~isempty(roiIndices)
            importFromVbFretDA(roiIndices);
        end
    end

    function exportAllFlaggedRoisForVbFretDA(obj,event)
        roiIndices = find(vertcat(rois.isFlagged) > 0);
        if ~isempty(roiIndices)
            exportForVbFretDA(roiIndices);
        end
    end

    function importIdealForAllFlaggedRoisFromVbFretDA(obj,event)
        roiIndices = find(vertcat(rois.isFlagged) > 0);
        if ~isempty(roiIndices)
            importFromVbFretDA(roiIndices);
        end
    end

    function exportForVbFretDA(roiIndices)
        disp('Exporting DA for vbFRET...');
        % Get save file.
        defaultSaveFilePath = '_vbFRET.mat';
        if ~isempty(strfind(roisFilePath,'.mat'))
            defaultSaveFilePath = [roisFilePath(1:end-4) '_vbFRET.mat'];
        end
        [file,path] = uiputfile(defaultSaveFilePath,'Export for vbFRET (*.mat).');
        if isequal(file,0)
            return
        end
        % Gen data array for vbFRET.
        data = {};
        n = 1;
        for k = 1:length(roiIndices)
            roiIndex = roiIndices(k);
            roi = rois(roiIndex);
            yDD = roi.zprojDD;
            yDA = roi.zprojDA;
            sDA = roi.splineDA;
            % Subtract donor cross talk from FRET channel.
            xtalk = str2num(get(crossTalkFractionEdit,'string'));
            if xtalk > 0
                yDA = yDA-(yDD.*xtalk);
            end
%             % Subtract non-constant baseline from FRET channel.
%             amp = roi.baselineDecayAmpDA;
%             if amp ~= 0
%                 tau = roi.baselineDecayTauDA;
%                 yDA = subtractBaselineDecay(yDA,amp,tau);
%             end
            % Subtract DA spline.
            if ~isempty(sDA) && length(sDA) == length(yDA)
                yDA = yDA-sDA;
            end
            % Chop trace a bit after bleaching.
            if ~isempty(roi.bleachFramesAA)
                bleachFrame = roi.bleachFramesAA(1);
                if bleachFrame+100 < length(yDA)
                    yDA = yDA(1:bleachFrame+100);
                end
            end
            % Compute artificial FRET to recapitulate DA channel.
            % FRET = A/(A+D)
            % Thus, if D = -A, then FRET is scaled version of A.
            acceptor = yDA-min(yDA);
            donor = -acceptor+max(acceptor);
            % Add to vbFRET array.
            data{n} = zeros(length(acceptor),2);
            data{n}(:,1) = donor;
            data{n}(:,2) = acceptor;
            n = n+1;
        end
        % Export to file.
        save([path file],'data');
        disp('... Finished exporting DA for vbFRET.');
    end

    function importFromVbFretDA(roiIndices)
        disp('Importing DA from vbFRET...');
        % Get open file.
        defaultOpenFilePath = '_vbFRET.mat';
        if ~isempty(strfind(roisFilePath,'.mat'))
            defaultOpenFilePath = [roisFilePath(1:end-4) '_vbFRET.mat'];
        end
        [file,path] = uigetfile(defaultOpenFilePath,'Import from vbFRET (*.mat).');
        if isequal(file,0)
            return
        end
        VB = load([path file]);
        for k = 1:length(roiIndices)
            ideal = VB.path{k};
            roiIndex = roiIndices(k);
            roi = rois(roiIndex);
            yDD = roi.zprojDD;
            yDA = roi.zprojDA;
            sDA = roi.splineDA;
            % Subtract donor cross talk from FRET channel.
            xtalk = str2num(get(crossTalkFractionEdit,'string'));
            if xtalk > 0
                yDA = yDA-(yDD.*xtalk);
            end
%             % Subtract non-constant baseline from FRET channel.
%             amp = roi.baselineDecayAmpDA;
%             if amp ~= 0
%                 tau = roi.baselineDecayTauDA;
%                 yDA = subtractBaselineDecay(yDA,amp,tau);
%             end
            % Subtract DA spline.
            if ~isempty(sDA) && length(sDA) == length(yDA)
                yDA = yDA-sDA;
            end
            % Chop trace a bit after bleaching.
            if ~isempty(roi.bleachFramesAA)
                bleachFrame = roi.bleachFramesAA(1);
                if bleachFrame+100 < length(yDA)
                    yDA = yDA(1:bleachFrame+100);
                end
            end
            acceptor = yDA-min(yDA);
            donor = -acceptor+max(acceptor);
            ideal = ideal.*(acceptor+donor);
            ideal = ideal+min(yDA);
            rois(roiIndex).idealDA = zeros(size(roi.zprojDA));
            rois(roiIndex).idealDA(1:length(ideal)) = ideal;
            if length(roi.zprojDA) > length(ideal)
                rois(roiIndex).idealDA(length(ideal)+1:end) = ideal(end);
            end
        end
        disp('... Finished importing DA from vbFRET.');
        goToRoi(); % Redraw.
    end

end
