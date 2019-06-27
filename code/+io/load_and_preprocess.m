function [stim, spks, params] = load_and_preprocess(fname, varargin)
% [stim, spks, params] = load_and_preprocess(fname, varargin)
% 

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

win = [0 0.3];
figure(1); clf
for iClust = 1:sp.nClusters
    spbinned = pdsa.binSpTimes(sp.spikeTimesSecs(sp.spikeClusters==sp.clusterId(iClust)), stimOnset, win, 0.01);

    subplot(4,sp.nClusters, (0:2)*sp.nClusters + iClust)
    imagesc(spbinned)
    subplot(4,sp.nClusters, 3*sp.nClusters + iClust)
    plot(nanmean(spbinned)) % used this to guess the optimal time lags
    
    spks(:,iClust) = nansum(spbinned(:,8),2); % guess the optimal time lags
end





    
    
    