
dataDir = uigetdir();

sessionList = dir(fullfile(dataDir, '*.mat'));

%%
iSess = 2;
[stim, spks, params] = io.load_and_preprocess(fullfile(dataDir, sessionList(iSess).name));

%% 
figure(1); clf
plot(spks)

%%
nStim = size(stim,2);
nTargets = unique(params.nFlashed);
[numShown, condId] = histc(params.nFlashed, [0 2 4 inf]);

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
        X = [stim(trialIx,:) ones(sum(trialIx),1)];
        sta = X'*R;
        ste = (X'*X + .1*eye(nStim+1))\(sta);
        sta(end) = [];
        
        
        figure(3); 
        
        subplot(1,nUnits, iUnit)
        plot(sta); hold on
        
        figure(2)
        subplot(nUnits, nConds, (iUnit-1)*nConds + iCond)
%         cmap = flipud(pdsa.cbrewer('seq', 'YlGnBu', 100));
        colormap(parula)
        staFull = match*sta;
        staFull(notshown) = nan;
        staFull = reshape(staFull, size(xx));
        plotWeights([xx(~notshown) yy(~notshown)], staFull(~notshown), 15*(rho(~notshown).^.55)); 
        
    end
end



%% Try fitting
% initlize with linear regression (across all conditions)
kUnit = 2;
R = spks(:,kUnit);
X = stim;

w0 = (X'*X)\(X'*R);

% fit the gain model

% we have to use expansive nonlinearities (exponential, powerlaw) sparingly
% because they can explode and get us into bad numerical territory

% specify the nonlinearities we'll fit with
f = @nlfuns.logexp1;
g = {@nlfuns.logexp1, @nlfuns.logistic};

% regularizaiton (this imposes smoothness in 2D space)
Cinv = eye(size(X,2)); %qfsmooth2nd(nx, ny); % 2nd derivative matrix
% Cinv = Cinv + eye(nx*ny); % Identity matrix --> ridge regression

% start by fitting the LNP model
mstruct = struct();
mstruct.neglogli = @regression.neglogli_poiss; % set the likelihood here
mstruct.logprior = @regression.logprior_Cinv;  % set the prior here
Xd = [X ones(size(X,1),1)]; % augment stimulus with a column of ones -- this creates a bias term

binSize = 1;
mstruct.liargs = {Xd, R, f, binSize};
mstruct.priargs = {blkdiag(Cinv, 0)};
fun = @(w) regression.neglogpost_GLM(w, 10, mstruct);

% options for fitting (turn on gradient checking)
% opts = optimset('GradObj', 'On', 'Hessian', 'Off', 'Display', 'none', 'MaxIter', 10e10);
% opts = optimoptions('fminunc','Algorithm','quasi-newton','Display', 'Final', 'SpecifyObjectiveGradient',true, 'MaxIter', 10e10);
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display', 'iter', 'SpecifyObjectiveGradient',true, 'MaxIter', 10e10, 'CheckGradients', true, 'OptimalityTolerance', 1e-9, 'FiniteDifferenceType', 'central');
w1 = fminunc(fun, [randn(size(w0)); mean(R)], opts);

figure(1); clf
subplot(1,3,1)
plotWeights([xx(~notshown) yy(~notshown)], w0, 15*rho(~notshown).^.5);
% plot(w0)
% imagesc(reshape(w0, [ny nx]))
title('Spike Triggred Average')

subplot(1,3,2)
plotWeights([xx(~notshown) yy(~notshown)], w1(1:end-1), 15*rho(~notshown).^.5);
% imagesc(reshape(w1(1:end-1), [ny nx]))
title('LNP')



figure(2); clf
plot(R, 'k');
hold on
RhatLN = f(Xd*w1)*binSize;
plot(RhatLN, 'r')
xlabel('Sample #')
ylabel('Spike Count')
legend('Data', 'LNP')

title(rsquared(R, RhatLN))

% initialize the gain field by fitting a hacky modulated poisson
Xd = [RhatLN(:) bsxfun(@times, X, RhatLN(:)) ones(numel(R),1)];

mstruct.liargs = {Xd, R, f, binSize};
mstruct.priargs = {blkdiag(.1, Cinv, 0)};

fun = @(w) regression.neglogpost_GLM(w, 1, mstruct);
wmod = fminunc(fun, [1; w1], opts);

figure(1)
subplot(1,3,3)
plotWeights([xx(~notshown) yy(~notshown)], wmod(2:end-1), 15*rho(~notshown).^.5);
title('Gain field Init')

% initialize filters
wGainInit =  wmod(2:end-1);
wRFInit = w1;

figure(2)
plot(f(Xd*wmod), 'g')
title([rsquared(R, RhatLN) rsquared(R, f(Xd*wmod))])
%% Fit gain model
% change the likelihood function to the modulated poisson model
mstruct.neglogli = @regression.neglogli_modulated_poiss;

% likelihood arguments: 
% {stimulus, spikes, subunit nonlinearity, spike, nonlinearity, bin size}
nTrials = size(X,1);
Xd = [X ones(nTrials,1)]; % back to the original stimulus
inds = 1:nTrials; % fit the RF using the one-target case
mstruct.liargs = {Xd(inds,:), R(inds), g{1}, f, binSize};
mstruct.priargs = {blkdiag(Cinv, 0.1)};

maxIter = 100; % maximum iterations to try

lambda = 1;

% set up the cost function
fun = @(w) regression.neglogpost_GLM(w, lambda, mstruct);

% wInitRF = w1; 
wInitRF = max(w1,0); % initialize with LNP, but only take positive values

% trying other initialization options
% wInitRF = randn(size(wInitRF));
% wHatRF = [Ke(:); rfOffset];
[wHatRF, fInit0] = fminunc(fun, wInitRF, opts);

figure(1)
subplot(2,3,1)
plotWeights([xx(~notshown) yy(~notshown)], wHatRF(1:end-1), 15*rho(~notshown).^.5);
title('RF')

% f(g{1}(

%%
% initialize denominator

% then generate the output of the RF, that gets passed into the likelihood
% function for the modulated poisson model
numerator = g{1}(Xd*wHatRF); % evaluate the signal in the numerator
% numerator = f(Xd*w1); % use the LNP model as the modulator signal

% swap out likelihood arguments for the divisive term 
mstruct.liargs = {Xd, R, g{2}, f, binSize, numerator};
fun = @(w) regression.neglogpost_GLM(w, lambda, mstruct);

% fun = @(w) regression.neglogli_modulated_poiss(w, X, R, g{2}, f, binSize, numerator);
% initialize with RF estimate? use STC?
% wInitGain = [Ks(:); dfOffset]; % if initializing with true filter
wInitGain = [wGainInit(:); 0.1]; %rand(size(wHatRF));
% wInitGain = rand(numel(wHatRF),1); % random initialization

[wHatGain, fInit] = fminunc(fun, wInitGain, opts);

figure(1)
subplot(2,3,1)
plotWeights([xx(~notshown) yy(~notshown)], wHatRF(1:end-1), 15*rho(~notshown).^.5);
title('RF')

subplot(2,3,4)
plotWeights([xx(~notshown) yy(~notshown)], wHatGain(1:end-1), 15*rho(~notshown).^.5);
title('div RF')
drawnow

% turn off gradient checking -- it fails sometimes  (I believe there are
% numerical errors from numerical chain rule the way I implemented it). If
% this continues
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display', 'None', 'SpecifyObjectiveGradient',true, 'MaxIter', 10e10, 'OptimalityTolerance', 1e-9);

modulator = g{2}(Xd*wHatGain);

figure(2); 
plot(f(numerator.*modulator)*binSize, 'b')
figure(3); clf
plot(numerator); hold on
plot(modulator)
%%
% initialize with our previous best estimates
wInitRF   = wHatRF;
wInitGain = wHatGain;

% random initializations -- probably won't work
% wInitRF   = randn(size(wHatRF));
% wInitGain = randn(size(wHatGain));

% fit iteratively
f1 = fInit;
figure(3); clf
subplot(1,2,1)
plot(wInitRF); hold on
subplot(1,2,2)
plot(wInitGain); hold on
drawnow
iter =0;
tol = 1e-3;

while iter < maxIter
    
    f0 = f1;
    modulator = g{2}(Xd*wInitGain);
    
    % fit numerator
    mstruct.liargs = {Xd, R, g{1}, f, binSize, modulator};
    mstruct.priargs = {blkdiag(Cinv, 1e10)};
    fun = @(w) regression.neglogpost_GLM(w, lambda, mstruct);
    
    [wHatRF,~] = fminunc(fun, wInitRF, opts);
    
    subplot(1,2,1)
    plot(wHatRF);
    title('RF')
    
    % fit denominator
    numerator = g{1}(Xd*wHatRF);
    mstruct.liargs = {Xd, R, g{2}, f, binSize, numerator};
    fun = @(w) regression.neglogpost_GLM(w, lambda, mstruct);

    [wHatGain,f1] = fminunc(fun, wInitGain, opts);
    
    % if the model stops improving, stop fitting
    if (f0 - f1) < tol 
        break
    end
    
    wInitGain = wHatGain;
    wInitRF = wHatRF;
    
    subplot(1,2,2)
    plot(wInitGain);
    title('div RF')
    drawnow
    
    
    iter = iter + 1; 
    fprintf('iter: %d, LL: %02.4f\n', iter, f1-f0)
end

% plot the outcome
wHatGain = wInitGain;
wHatRF = wInitRF;

figure(4); clf

subplot(2,2,3)
plotWeights([xx(~notshown) yy(~notshown)], wHatRF(1:end-1), 15*rho(~notshown).^.5);
title('K_e fit')

subplot(2,2,4)
plotWeights([xx(~notshown) yy(~notshown)], wHatGain(1:end-1), 15*rho(~notshown).^.5);
title('K_s fit')

set(gcf, 'PaperSize', [5 5], 'PaperPosition', [0 0 5 5])
% saveas(gcf, fullfile(figDir, 'divSfits.png'))

Rhat = f(g{1}(Xd*wHatRF).*g{2}(Xd*wHatGain))*binSize;

figure(2); clf
plot(R, 'k'); hold on
plot(Rhat, 'g')
xlabel('Sample #')
ylabel('Spike Count')
legend('Data', 'divS')
title(rsquared(R, Rhat))


set(gcf, 'PaperSize', [6 3], 'PaperPosition', [0 0 6 3])
% saveas(gcf, fullfile(figDir, 'divSSTA.png'))
