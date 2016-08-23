function [analysis,rois] = analyzeRois(rois,visualize)
%%
% Author: Marcel Goldschen-Ohm
% Email: marcel.goldschen@gmail.com
analysis = struct();

%% Unique ligand values.
fcAMP_uM = unique(vertcat(rois.fcAMP_uM));
analysis.fcAMP_uM = fcAMP_uM;

%% Find events.
for k = 1:length(rois)
    roi = rois(k);
    events = [];
    bf = roi.bleachFramesAA(1);
    if bf > 3
        ideal = roi.idealDA(1:bf-2);
        events = [events;[1,ideal(1)]];
        for j = 2:length(ideal)
            if ideal(j) ~= events(end,2)
                events = [events;[j,ideal(j)]];
            end
        end
        events = [events;[bf-1,-1]];
    end
    rois(k).events = events;
end
clear k j roi events bf ideal

%% Find dwell times and bound probability vs [fcAMP].
analysis.unboundTimes_s = {};
analysis.boundTimes_s = {};
analysis.unboundBoundTimePairs_s = {};
analysis.boundUnboundTimePairs_s = {};
analysis.unboundUnboundTimePairs_s = {};
analysis.boundBoundTimePairs_s = {};
for j = 1:length(fcAMP_uM)
    analysis.unboundTimes_s{j} = [];
    analysis.boundTimes_s{j} = [];
    analysis.unboundBoundTimePairs_s{j} = [];
    analysis.boundUnboundTimePairs_s{j} = [];
    analysis.unboundUnboundTimePairs_s{j} = [];
    analysis.boundBoundTimePairs_s{j} = [];
end
for k = 1:length(rois)
    roi = rois(k);
    % Find fcAMP index.
    j = find(fcAMP_uM == roi.fcAMP_uM);
    % Idealization.
    bleachFrame = roi.bleachFramesAA(1);
    ideal = roi.idealDA(1:bleachFrame-1);
    % Bound probability and stuff.
    rois(k).boundProbability = double(sum(ideal))/length(ideal);
    rois(k).totalUnboundTime_s = double(length(find(ideal == 0)))*0.2;
    rois(k).totalBoundTime_s = double(length(find(ideal > 0)))*0.2;
    % Dwell times for this ROI.
    rois(k).avgUnboundTime_s = 0;
    rois(k).avgBoundTime_s = 0;
    if size(roi.events,1) > 2
        frames = diff(roi.events(:,1));
        states = roi.events(1:end-1,2);
        if length(states) > 2
            % Drop first and last event.
            frames = frames(2:end-1);
            states = states(2:end-1);
            % Bleach frame.
            unbound = find(states == 0);
            bound = find(states == 1);
            if ~isempty(unbound)
                analysis.unboundTimes_s{j} = [analysis.unboundTimes_s{j};double(frames(unbound)).*0.2];
                rois(k).avgUnboundTime_s = mean(double(frames(unbound)).*0.2);
            end
            if ~isempty(bound)
                analysis.boundTimes_s{j} = [analysis.boundTimes_s{j};double(frames(bound)).*0.2];
                rois(k).avgBoundTime_s = mean(double(frames(bound)).*0.2);
            end
            numDwells = length(frames);
            % 1st order correlations.
            if numDwells >= 2
                numPairs = (numDwells-mod(numDwells,2))/2;
                pairs = zeros(numPairs,2);
                pairs(:,1) = frames(1:2:2*numPairs);
                pairs(:,2) = frames(2:2:2*numPairs);
                if states(1) == 0
                    analysis.unboundBoundTimePairs_s{j} = [analysis.unboundBoundTimePairs_s{j};double(pairs).*0.2];
                else
                    analysis.boundUnboundTimePairs_s{j} = [analysis.boundUnboundTimePairs_s{j};double(pairs).*0.2];
                end
                if numDwells >= 3
                    numPairs = (numDwells-1-mod(numDwells-1,2))/2;
                    pairs = zeros(numPairs,2);
                    pairs(:,1) = frames(2:2:1+2*numPairs);
                    pairs(:,2) = frames(3:2:1+2*numPairs);
                    if states(2) == 0
                        analysis.unboundBoundTimePairs_s{j} = [analysis.unboundBoundTimePairs_s{j};double(pairs).*0.2];
                    else
                        analysis.boundUnboundTimePairs_s{j} = [analysis.boundUnboundTimePairs_s{j};double(pairs).*0.2];
                    end
                end
            end
            % 2nd order correlations.
            if numDwells >= 3
                numPairs = int32(floor(double(numDwells-1)/2));
                pairs = zeros(numPairs,2);
                pairs(:,1) = frames(1:2:2*numPairs);
                pairs(:,2) = frames(3:2:1+2*numPairs);
                if states(1) == 0
                    analysis.unboundUnboundTimePairs_s{j} = [analysis.unboundUnboundTimePairs_s{j};double(pairs).*0.2];
                else
                    analysis.boundBoundTimePairs_s{j} = [analysis.boundBoundTimePairs_s{j};double(pairs).*0.2];
                end
                if numDwells >= 4
                    numPairs = int32(floor(double(numDwells-2)/2));
                    pairs = zeros(numPairs,2);
                    pairs(:,1) = frames(2:2:1+2*numPairs);
                    pairs(:,2) = frames(4:2:2+2*numPairs);
                    if states(2) == 0
                        analysis.unboundUnboundTimePairs_s{j} = [analysis.unboundUnboundTimePairs_s{j};double(pairs).*0.2];
                    else
                        analysis.boundBoundTimePairs_s{j} = [analysis.boundBoundTimePairs_s{j};double(pairs).*0.2];
                    end
                end
            end
        end
    end
end
clear k j roi bleachFrame bound unbound frames states numDwells numPairs pairs

%% Compute bound probability vs [fcAMP].
analysis.boundProbability = zeros(size(fcAMP_uM));
analysis.boundProbabilitySEM = zeros(size(fcAMP_uM));
analysis.nullProbability = zeros(size(fcAMP_uM));
% Pnull = [.6789,.4196,.1875,.0558,.0080];
for j = 1:length(fcAMP_uM)
    idx = find(vertcat(rois.fcAMP_uM) == fcAMP_uM(j));
    totalUnboundTime_s = sum(vertcat(rois(idx).totalUnboundTime_s));
    totalBoundTime_s = sum(vertcat(rois(idx).totalBoundTime_s));
    analysis.boundProbability(j) = totalBoundTime_s/(totalBoundTime_s+totalUnboundTime_s);
    % analysis.boundProbability(j) = totalBoundTime_s/(totalBoundTime_s+totalUnboundTime_s*(1+Pnull(j)));
    perRoiBoundProb = vertcat(rois(idx).boundProbability);
    perRoiBoundProb = perRoiBoundProb(find(~isnan(perRoiBoundProb)));
    if ~isempty(perRoiBoundProb)
        analysis.boundProbabilitySEM(j) = std(perRoiBoundProb) / sqrt(length(perRoiBoundProb));
    end
    analysis.nullProbability(j) = double(length(find(perRoiBoundProb == 0)))/length(perRoiBoundProb);
end
clear j idx totalUnboundTime_s totalBoundTime_s perRoiBoundProb

%% Fit binding curve Bmax * X / (X + Kd).
analysis.bindingCurve = @(p,x) (p(1).*x)./(x+p(2));
costFunc = @(p) sum((analysis.bindingCurve(p,fcAMP_uM) - analysis.boundProbability).^2 ./ analysis.boundProbabilitySEM);
Bmax = 1;
Kd = 1;
p = fminsearch(costFunc,[Bmax,Kd]);
analysis.Bmax = p(1);
analysis.Kd = p(2);
clear p Bmax Kd costFunc

%% Plot bound probability vs [fcAMP] with fitted binding curve.
if visualize
    figure; hold on;
    errorbar(fcAMP_uM,analysis.boundProbability,analysis.boundProbabilitySEM,'bo');
    plot(0:0.1:10,analysis.bindingCurve([analysis.Bmax,analysis.Kd],0:0.1:10),'r-');
    xlabel('[fcAMP] (uM)');
    ylabel('Bound Probability');
    title(['Bmax = ' num2str(analysis.Bmax) ', Kd = ' num2str(analysis.Kd) ' uM']);
end

%% Bleach times.
analysis.bleachTimes_s = zeros(size(rois));
for k = 1:length(rois)
    analysis.bleachTimes_s(k) = double(rois(k).bleachFramesAA(1))*0.2;
end
% Fit with exponential.
[phat,pci] = mle(analysis.bleachTimes_s,'distribution','exp');
analysis.bleachMonoExpTau_s = phat;
analysis.bleachMonoExpTauConfidenceInterval_s = pci;
% Fit with stretched exponential.
stretchedExpPdf = @(t,tau,beta) 1/(tau*gamma(1+1/beta)).*exp(-(t./tau).^beta);
[phat,pci] = mle(analysis.bleachTimes_s,'pdf',stretchedExpPdf,'start',[analysis.bleachMonoExpTau_s,0.5],'lowerbound',[0.2,0],'upperbound',[1000,1]);
analysis.bleachStretchedExpTau_s = phat(1);
analysis.bleachStretchedExpBeta = phat(2);
analysis.bleachStretchedExpTauConfidenceInterval_s = pci(:,1);
analysis.bleachStretchedExpBetaConfidenceInterval = pci(:,2);
clear k phat pci

%% Plot bleach times.
if visualize
    figure; hold on;
    bins = linspace(0,200,50);
    counts = hist(analysis.bleachTimes_s,bins);
    x = bins(1:end-1);
    y = counts(1:end-1);
    stairs(x,y,'b-');
    plot(x,exppdf(x,analysis.bleachMonoExpTau_s).*trapz(x,y),'r-');
    plot(x,stretchedExpPdf(x,analysis.bleachStretchedExpTau_s,analysis.bleachStretchedExpBeta).*trapz(x,y),'r--');
    ylabel('Counts');
    xlabel('Bleach Time (s)');
    legend({...
        'Bleach Times'...
        ,['tau = ' num2str(analysis.bleachMonoExpTau_s,3) ' sec']...
        ,['tau = ' num2str(analysis.bleachStretchedExpTau_s,3) ' sec, beta = ' num2str(analysis.bleachStretchedExpBeta,2)]...
        });
    clear bins counts x y
end

%% Fit dwell time distributions with monoexponentials.
analysis.unboundMonoExpTaus_s = zeros(length(fcAMP_uM),1);
analysis.boundMonoExpTaus_s = zeros(length(fcAMP_uM),1);
analysis.unboundMonoExpTauConfidenceIntervals_s = zeros(length(fcAMP_uM),2);
analysis.boundMonoExpTauConfidenceIntervals_s = zeros(length(fcAMP_uM),2);
for j = 1:length(fcAMP_uM)
    [phat,pci] = mle(analysis.unboundTimes_s{j},'distribution','exp');
    analysis.unboundMonoExpTaus_s(j) = phat;
    analysis.unboundMonoExpTauConfidenceIntervals_s(j,:) = pci(:,1);
    [phat,pci] = mle(analysis.boundTimes_s{j},'distribution','exp');
    analysis.boundMonoExpTaus_s(j) = phat;
    analysis.boundMonoExpTauConfidenceIntervals_s(j,:) = pci(:,1);
end
clear j phat pci

%% Fit dwell time distributions with stretched exponentials.
stretchedExpPdf = @(t,tau,beta) 1/(tau*gamma(1+1/beta)).*exp(-(t./tau).^beta);
analysis.unboundStretchedExpTaus_s = zeros(length(fcAMP_uM),1);
analysis.unboundStretchedExpBetas = zeros(length(fcAMP_uM),1);
analysis.boundStretchedExpTaus_s = zeros(length(fcAMP_uM),1);
analysis.boundStretchedExpBetas = zeros(length(fcAMP_uM),1);
analysis.unboundStretchedExpTauConfidenceIntervals_s = zeros(length(fcAMP_uM),2);
analysis.unboundStretchedExpBetaConfidenceIntervals = zeros(length(fcAMP_uM),2);
analysis.boundStretchedExpTauConfidenceIntervals_s = zeros(length(fcAMP_uM),2);
analysis.boundStretchedExpBetaConfidenceIntervals = zeros(length(fcAMP_uM),2);
for j = 1:length(fcAMP_uM)
    [phat,pci] = mle(analysis.unboundTimes_s{j},'pdf',stretchedExpPdf,'start',[analysis.unboundMonoExpTaus_s(j),0.5],'lowerbound',[0.01,0],'upperbound',[1000,1]);
    analysis.unboundStretchedExpTaus_s(j) = phat(1);
    analysis.unboundStretchedExpBetas(j) = phat(2);
    analysis.unboundStretchedExpTauConfidenceIntervals_s(j,:) = pci(:,1);
    analysis.unboundStretchedExpBetaConfidenceIntervals(j,:) = pci(:,2);
    [phat,pci] = mle(analysis.boundTimes_s{j},'pdf',stretchedExpPdf,'start',[analysis.boundMonoExpTaus_s(j),0.5],'lowerbound',[0.01,0],'upperbound',[1000,1]);
    analysis.boundStretchedExpTaus_s(j) = phat(1);
    analysis.boundStretchedExpBetas(j) = phat(2);
    analysis.boundStretchedExpTauConfidenceIntervals_s(j,:) = pci(:,1);
    analysis.boundStretchedExpBetaConfidenceIntervals(j,:) = pci(:,2);
end
clear j phat pci

%% Fit dwell time distributions with with biexponentials.
twoExpPdf = @(t,A1,tau1,tau2) A1.*exppdf(t,tau1) + (1-A1).*exppdf(t,tau2);
analysis.unboundBiExpAmps = zeros(length(fcAMP_uM),2);
analysis.unboundBiExpTaus_s = zeros(length(fcAMP_uM),2);
analysis.unboundBiExpAmpConfidenceIntervals = zeros(length(fcAMP_uM),4);
analysis.unboundBiExpTauConfidenceIntervals_s = zeros(length(fcAMP_uM),4);
analysis.boundBiExpAmps = zeros(length(fcAMP_uM),2);
analysis.boundBiExpTaus_s = zeros(length(fcAMP_uM),2);
analysis.boundBiExpAmpConfidenceIntervals = zeros(length(fcAMP_uM),4);
analysis.boundBiExpTauConfidenceIntervals_s = zeros(length(fcAMP_uM),4);
for j = 1:length(fcAMP_uM)
    [phat,pci] = mle(analysis.unboundTimes_s{j},'pdf',twoExpPdf,'start',[0.5,analysis.unboundStretchedExpTaus_s(j),analysis.unboundMonoExpTaus_s(j)],'lowerbound',[0,0.01,0.2],'upperbound',[1,1000,1000]);
    analysis.unboundBiExpAmps(j,:) = [phat(1),1-phat(1)];
    analysis.unboundBiExpTaus_s(j,:) = [phat(2),phat(3)];
    analysis.unboundBiExpAmpConfidenceIntervals(j,1:2) = pci(:,1);
    analysis.unboundBiExpAmpConfidenceIntervals(j,3:4) = [1-pci(2,1),1-pci(1,1)];
    analysis.unboundBiExpTauConfidenceIntervals_s(j,1:2) = pci(:,2);
    analysis.unboundBiExpTauConfidenceIntervals_s(j,3:4) = pci(:,3);
    [phat,pci] = mle(analysis.boundTimes_s{j},'pdf',twoExpPdf,'start',[0.5,analysis.boundStretchedExpTaus_s(j),analysis.boundMonoExpTaus_s(j)],'lowerbound',[0,0.01,0.2],'upperbound',[1,1000,1000]);
    analysis.boundBiExpAmps(j,:) = [phat(1),1-phat(1)];
    analysis.boundBiExpTaus_s(j,:) = [phat(2),phat(3)];
    analysis.boundBiExpAmpConfidenceIntervals(j,1:2) = pci(:,1);
    analysis.boundBiExpAmpConfidenceIntervals(j,3:4) = [1-pci(2,1),1-pci(1,1)];
    analysis.boundBiExpTauConfidenceIntervals_s(j,1:2) = pci(:,2);
    analysis.boundBiExpTauConfidenceIntervals_s(j,3:4) = pci(:,3);
end
clear j phat pci

%% Sort biexponential time constants into major and minor unbound amplitudes and fast and slow bound times.
for j = 1:length(fcAMP_uM)
    if analysis.unboundBiExpAmps(j,1) < analysis.unboundBiExpAmps(j,2)
        analysis.unboundBiExpAmps(j,:) = [analysis.unboundBiExpAmps(j,2),analysis.unboundBiExpAmps(j,1)];
        analysis.unboundBiExpTaus_s(j,:) = [analysis.unboundBiExpTaus_s(j,2),analysis.unboundBiExpTaus_s(j,1)];
        analysis.unboundBiExpAmpConfidenceIntervals(j,:) = ...
            [analysis.unboundBiExpAmpConfidenceIntervals(j,3:4),analysis.unboundBiExpAmpConfidenceIntervals(j,1:2)];
        analysis.unboundBiExpTauConfidenceIntervals_s(j,:) = ...
            [analysis.unboundBiExpTauConfidenceIntervals_s(j,3:4),analysis.unboundBiExpTauConfidenceIntervals_s(j,1:2)];
    end
    if analysis.boundBiExpTaus_s(j,1) > analysis.boundBiExpTaus_s(j,2)
        analysis.boundBiExpAmps(j,:) = [analysis.boundBiExpAmps(j,2),analysis.boundBiExpAmps(j,1)];
        analysis.boundBiExpTaus_s(j,:) = [analysis.boundBiExpTaus_s(j,2),analysis.boundBiExpTaus_s(j,1)];
        analysis.boundBiExpAmpConfidenceIntervals(j,:) = ...
            [analysis.boundBiExpAmpConfidenceIntervals(j,3:4),analysis.boundBiExpAmpConfidenceIntervals(j,1:2)];
        analysis.boundBiExpTauConfidenceIntervals_s(j,:) = ...
            [analysis.boundBiExpTauConfidenceIntervals_s(j,3:4),analysis.boundBiExpTauConfidenceIntervals_s(j,1:2)];
    end
end
clear j

%% Plot dwell times with fits.
if visualize
    figure;
    for j = 1:length(fcAMP_uM)
        subplot(length(fcAMP_uM),2,(j-1)*2+1); cla; hold on;
        bins = linspace(0,30,75);
        counts = hist(analysis.unboundTimes_s{j},bins);
        counts = counts./trapz(bins,counts);
        x = bins(1:end-1);
        y = counts(1:end-1);
        stairs(x,y,'b-')
        plot(x,exppdf(x,analysis.unboundMonoExpTaus_s(j)).*trapz(x,y),'-','color',[.5,.5,.5]);
        plot(x,twoExpPdf(x,analysis.unboundBiExpAmps(j,1),analysis.unboundBiExpTaus_s(j,1),analysis.unboundBiExpTaus_s(j,2)).*trapz(x,y),'r-');
        xlabel('Unbound (s)');
        ylabel({[num2str(fcAMP_uM(j)) ' uM'],'Counts'});

        subplot(length(fcAMP_uM),2,(j-1)*2+2); cla; hold on;
        bins = linspace(0,10,50);
        counts = hist(analysis.boundTimes_s{j},bins);
        counts = counts./trapz(bins,counts);
        x = bins(1:end-1);
        y = counts(1:end-1);
        stairs(x,y,'b-');
        plot(x,exppdf(x,analysis.boundMonoExpTaus_s(j)).*trapz(x,y),'-','color',[.5,.5,.5]);
        plot(x,twoExpPdf(x,analysis.boundBiExpAmps(j,1),analysis.boundBiExpTaus_s(j,1),analysis.boundBiExpTaus_s(j,2)).*trapz(x,y),'r-');
        xlabel('Bound (s)');
        ylabel({[num2str(fcAMP_uM(j)) ' uM'],'Counts'}); 
    end
end

%% Plot summaries of dwell time biexponential fits.
if visualize
    figure;
    subplot(2,2,1); cla; hold on;
    errorbar(fcAMP_uM,analysis.unboundBiExpAmps(:,2),analysis.unboundBiExpAmps(:,2)-analysis.unboundBiExpAmpConfidenceIntervals(:,3),analysis.unboundBiExpAmpConfidenceIntervals(:,4)-analysis.unboundBiExpAmps(:,2),'ro-');
    errorbar(fcAMP_uM,analysis.unboundBiExpAmps(:,1),analysis.unboundBiExpAmps(:,1)-analysis.unboundBiExpAmpConfidenceIntervals(:,1),analysis.unboundBiExpAmpConfidenceIntervals(:,2)-analysis.unboundBiExpAmps(:,1),'bo-');
    xlim([0,10.1]);
    ylim([0,1]);
    xlabel('fcAMP (uM)');
    ylabel('Unbound Amp');

    subplot(2,2,2); cla; hold on;
    errorbar(fcAMP_uM,analysis.boundBiExpAmps(:,2),analysis.boundBiExpAmps(:,2)-analysis.boundBiExpAmpConfidenceIntervals(:,3),analysis.boundBiExpAmpConfidenceIntervals(:,4)-analysis.boundBiExpAmps(:,2),'ro-');
    errorbar(fcAMP_uM,analysis.boundBiExpAmps(:,1),analysis.boundBiExpAmps(:,1)-analysis.boundBiExpAmpConfidenceIntervals(:,1),analysis.boundBiExpAmpConfidenceIntervals(:,2)-analysis.boundBiExpAmps(:,1),'bo-');
    xlim([0,10.1]);
    ylim([0,1]);
    xlabel('fcAMP (uM)');
    ylabel('Bound Amp');

    subplot(2,2,3); cla; hold on;
    errorbar(fcAMP_uM,analysis.unboundBiExpTaus_s(:,2),analysis.unboundBiExpTaus_s(:,2)-analysis.unboundBiExpTauConfidenceIntervals_s(:,3),analysis.unboundBiExpTauConfidenceIntervals_s(:,4)-analysis.unboundBiExpTaus_s(:,2),'ro-');
    errorbar(fcAMP_uM,analysis.unboundBiExpTaus_s(:,1),analysis.unboundBiExpTaus_s(:,1)-analysis.unboundBiExpTauConfidenceIntervals_s(:,1),analysis.unboundBiExpTauConfidenceIntervals_s(:,2)-analysis.unboundBiExpTaus_s(:,1),'bo-');
    xlim([0,10.1]);
    xlabel('fcAMP (uM)');
    ylabel('Unbound Tau (s)');

    subplot(2,2,4); cla; hold on;
    errorbar(fcAMP_uM,analysis.boundBiExpTaus_s(:,2),analysis.boundBiExpTaus_s(:,2)-analysis.boundBiExpTauConfidenceIntervals_s(:,3),analysis.boundBiExpTauConfidenceIntervals_s(:,4)-analysis.boundBiExpTaus_s(:,2),'ro-');
    errorbar(fcAMP_uM,analysis.boundBiExpTaus_s(:,1),analysis.boundBiExpTaus_s(:,1)-analysis.boundBiExpTauConfidenceIntervals_s(:,1),analysis.boundBiExpTauConfidenceIntervals_s(:,2)-analysis.boundBiExpTaus_s(:,1),'bo-');
    xlim([0,10.1]);
    xlabel('fcAMP (uM)');
    ylabel('Bound Tau (s)');
end
