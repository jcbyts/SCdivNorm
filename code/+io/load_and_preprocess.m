function [stim, spks, params] = load_and_preprocess(fname, varargin)
% [stim, spks, params] = load_and_preprocess(fname, varargin)
% This is just a quick import function
% make a matrix for the stim and count spikes

if ~exist(fname, 'file')
    if ~isempty(getpref('scDivNorm'))
        dataDir = getpref('scDivNorm', 'dataDir');
    else
        dataDir = uigetdir();
    end
    fname = fullfile(dataDir, fname);
end


load(fname)

% From Leor:
% - Take trials that are:
% bb.mFlash.isGoodTrial==1
% (though in theory, you could also use the bad trials pre-fixation break, but lets not get greedy now)
% 
% - stimOn times are in:
% bb.mFlash.timing.stimOn
% which is a cellarray. Each trial is a cell. 
% The number of elements per cell corresponds to the number of presentations in that trial.
% 
% Then it gets stupid hairy....
% your stimulus location (in cartesian and polar) are in a struct called "flash", and you can index into it with trial index or presentation index:
% bb.mFlash.stim(iTrial).flash(iPresentaion)
% 
% That's basically it for behavior. 
% 
% Spikes:
% I made an sp-like structure (even though this is single channel, sorted by plxOfs), for your-sp-convenience. 
% sp.



flashCtr = 1;

nTrials = numel(bb.mFlash.timing.stimOn);
nFlashTotal = sum(cellfun(@numel, bb.mFlash.timing.stimOn));
stimOnset = nan(nFlashTotal,1);

stimTh = cell(nFlashTotal,1);
stimRho = cell(nFlashTotal,1);

for iTrial = 1:nTrials
    
    onsetTimes = bb.mFlash.timing.stimOn{iTrial};
    nPresentations = numel(onsetTimes);
    if ~isnan(bb.mFlash.timing.fixBreak(iTrial))
        nPresentations = find(onsetTimes < bb.mFlash.timing.fixBreak(iTrial), 1, 'last');
    end
    
    for iPresentation = 1:nPresentations
        if ischar(bb.mFlash.stim(iTrial).flash(iPresentation).theta)
            continue
        end
        stimTh{flashCtr} = bb.mFlash.stim(iTrial).flash(iPresentation).theta(:);
        stimRho{flashCtr} = bb.mFlash.stim(iTrial).flash(iPresentation).radius(:);
        stimOnset(flashCtr) = onsetTimes(iPresentation);
        flashCtr = flashCtr + 1;
    end
    
end
    

stimPolar = [cell2mat(stimTh), cell2mat(stimRho)];

stimPolar = unique(stimPolar(~any(isnan(stimPolar),2),:), 'rows');
nStim = size(stimPolar,1);

% bin up the stimulus
stim = zeros(flashCtr-1, nStim);
for iFlash = 1:flashCtr-1
    stim(iFlash,:) = sum(double(stimPolar(:,1)==stimTh{iFlash}' & stimPolar(:,2)==stimRho{iFlash}'),2);    
end

stimOnset = stimOnset(1:flashCtr-1);
params.nFlashed = sum(stim,2);
params.stimPolarCoords = stimPolar;
params.thetas = unique(params.stimPolarCoords(:,1));
params.rhos   = unique(params.stimPolarCoords(:,2));
params.nThetas = numel(params.thetas);
params.nRhos  = numel(params.rhos);

spks = zeros(flashCtr-1,sp.nClusters);

% This is all for testing out a spatiotemporal RF, but probably shouldn't
% be used
% %% build desgn matrix and check that things worked
% t0 = min(stimOnset);
% t1 = max(stimOnset) + 1;
% binSize = 0.02;
% timeBins = t0:binSize:t1;
% T = numel(timeBins);
% nkt = ceil(.3/binSize);
% 
% flashBinned = histc(stimOnset, timeBins);
% 
% spbinned = nan(T, sp.nClusters);
% for iClust = 1:sp.nClusters
%     spbinned(:,iClust) = histc(sp.spikeTimesSecs(sp.spikeClusters==sp.clusterId(iClust)), timeBins);
% end
% 
% X = zeros(T, size(stim,2));
% inds = find(flashBinned);
% for i = 1:numel(inds)
%     X(inds(i),:) = stim(i,:);
% end
% 
% Xd = makeStimRows(X, nkt);
% %%
% y = spbinned - mean(spbinned);
% sta = Xd'*y;
% 
% %%
% figure(1); clf

% ax = pdsa.tight_subplot(sx, sy, 0.01, 0.01);
% for iClust = 1:sp.nClusters
%     set(1, 'CurrentAxes', ax(iClust))
% %     imagesc(reshape(sta(:,iClust), nkt, []))
%     plot(mean(reshape(sta(:,iClust), nkt, []), 2)); hold on
%     plot(xlim, [0 0], 'k--')
%     axis off
% end
% 
% 
% 
% %%
% k_rank = 1;
% [what1map,wt1map,wx1map] = bilinearRegress_coordAscent_xtMAP(X,spbinned,num_lags,k_rank,Ctinv,Cxinv,opts);
% plot(spbinned)
% %%

%% count spikes in a window after the flash
win = [0 0.3];
figure(1); clf
sx = ceil(sqrt(sp.nClusters));
sy = round(sqrt(sp.nClusters));
ax = pdsa.tight_subplot(sx, sy, 0.02, 0.1);
binSize = 0.01;
countingWindow = [.04 .15]; % bins after onset

for iClust = 1:sp.nClusters
    [spbinned, bcenters] = pdsa.binSpTimes(sp.spikeTimesSecs(sp.spikeClusters==sp.clusterId(iClust)), stimOnset, win, binSize);

    set(gcf, 'CurrentAxes', ax(iClust))
    plot(bcenters, nanmean(spbinned)/binSize) % used this to guess the optimal time lags
    hold on
    if iClust <= (sy-1) * sx
        set(gca, 'XTickLabel', '')
    end
    fill(countingWindow([1 1 2 2]), [ylim, fliplr(ylim)], 'k', 'FaceAlpha', .2, 'EdgeColor', 'none' )
    spks(:,iClust) = nansum(spbinned(:,5:8),2); % guess the optimal time lags
    title(iClust)
end

pdsa.fixfigure(gcf, 8, [12 12])





    
    
    