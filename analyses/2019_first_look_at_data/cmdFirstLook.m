
dataDir = uigetdir();

sessionList = dir(fullfile(dataDir, '*.mat'));

%%
iSess = 1;
[stim, spks, params] = io.load_and_preprocess(fullfile(dataDir, sessionList(iSess).name));


%%
nTargets = unique(params.nFlashed);
[numShown, condId] = histc(params.nFlashed, [0 2 4 inf]);
numShown
% condId = ones(size(condId));
conds = unique(condId);
nConds = numel(conds);
nUnits = size(spks,2);


[th,rho] = meshgrid(params.thetas, params.rhos);

match = th(:)==params.stimPolarCoords(:,1)' & rho(:)==params.stimPolarCoords(:,2)';
notshown = sum(match,2)==0;

xx = rho.*cosd(th);
yy = rho.*sind(th);


% plot the responses, check that the STA changes with # targets
figure(2); clf
figure(3); clf
for iCond = 1:nConds
    
    trialIx = condId==conds(iCond);
    for iUnit = 1:nUnits
        
        R = spks(trialIx, iUnit);
        R = R - mean(R); % subtract baseline
        X = stim(trialIx,:);
        sta = X'*R;
        sta = sta./sum(X)';%         Xsum(X);
        sta(isnan(sta)) = 0;
        figure(3); 
        subplot(1,nUnits, iUnit)
        plot(sta); hold on
        figure(2)
        subplot(nUnits, nConds, (iUnit-1)*nConds + iCond)
        staFull = match*sta;
        staFull(notshown) = nan;
        h = pcolor(xx, yy, reshape(staFull, size(xx))); 
        h.LineStyle = 'none';
        h.AlphaData = reshape(staFull, size(xx));
%         axis 
    end
end

%%
n = 6;
r = (0:n)'/n;
theta = pi*(-n:n)/n;
X = r*cos(theta);
Y = r*sin(theta);
C = r*cos(2*theta);
pcolor(X,Y,C)
axis equal tight 