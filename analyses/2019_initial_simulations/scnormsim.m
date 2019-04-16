% Let's try to see if this makes any sense as an approach
% First off, we want to simulate some responses from a neuron that has
% divisive suppression. 

figDir = fullfile(fileparts(which('addSCdivNorm')), 'figures');

% Hyperparameters: play with them
% these parameters apply regularization -- currently they will smooth and
% shrink the weights. We will have to select these from data using cross
% validation. It looks like the model is pretty quick to fit, so that
% should be so bad
rho = 0.1; 
rho_sup = 2; % hyper parameter for the divisive RF
nSamples = 500; % number samples per condition
%% generate some data
nx = 8; % number of spatial bins in the x-dimension
ny = 10; % number of spatial bins in the y-dimension

% parameters to govern the shape of the receptive/divisive fields and the
% sensitivity of the neuron
rfCenter = [6 4];
dfCenter = [2 5];
rfSigma = 1; % standard deviation of RF
dfSigma = 2; % standard deviation of DF
rfScale = 40; % scales up RF drive
dfScale = -2; % scales divisive drive
rfOffset = 4;
dfOffset = 1.5; % shifts offset point of the divisive drive


% build the receptive and divisive fields 
xax = 1:nx;
yax = 1:ny;
[xx,yy] = meshgrid(xax, yax);

Ke = rfScale * exp(-.5 * ((xx - rfCenter(1)).^2 + (yy - rfCenter(2)).^2)/rfSigma^2);
Ks = dfScale*exp(-.5 * ((xx - dfCenter(1)).^2 + (yy - dfCenter(2)).^2)/dfSigma^2);

binSize    = .1; % 100 ms

% specify the nonlinearities
g_e     = @nlfuns.threshLinear;
g_s     = @nlfuns.logistic;
f       = @nlfuns.logexp1;

% Build the stimulus
nTargets = [1 2 4 8]; % number of targets in each condition
ndim = nx*ny;

figure(1); clf
subplot(1,2,1)
imagesc(Ke)
title('Receptive Field')

subplot(1,2,2)
imagesc(Ks)
title('Divisive Field')

% generate stimulus
nConds = numel(nTargets);
X = zeros(nSamples*nConds, ndim); % this is how we represent the targets on the stimulus
targetCond = zeros(nSamples*nConds,1);
for j = 1:nConds
    for i = 1:nSamples
        stimId = randi(ndim, 1, nTargets(j));
        X((j-1)*nSamples+i,stimId) = 1;
        targetCond((j-1)*nSamples+i) = nTargets(j);
    end
end

% generate spiking responses
rfDrive = X*Ke(:)+rfOffset;
dfDrive = X*Ks(:)+dfOffset;
lambda = f( g_e(rfDrive) .* g_s(dfDrive));

R = (lambda*binSize);
R = poissrnd(R); % comment out if you want to check the noiseless case

% plot the responses, check that the STA changes with # targets
figure(2); clf
for i = nTargets
    inds = find(targetCond==i);
    subplot(2,nConds,1:nConds)
    plot(inds, R(inds), '-'); hold on
    plot(inds([1 end]), mean(R(inds))*[1 1], 'k-')
    title('the "Data"')
    ylabel('Spike Count')
    xlabel('Sample #')
    
    subplot(2,nConds,nConds+find(i==nTargets))
    sta = simpleSTC(X(inds,:),R(inds), 1);
%     sta = (X(inds,:)'*R(inds))/sum(R(inds));
    sta = reshape(sta, [ny nx]);
    imagesc(sta)
    title(sprintf('STA (%d targs)', i))
    
end
set(gcf, 'PaperSize', [6 3], 'PaperPosition', [0 0 6 3])
saveas(gcf, fullfile(figDir, 'simulationSTA.png'))
%% Try fitting
% initlize with linear regression (across all conditions)
w0 = (X'*X)\(X'*R);




% fit the gain model

% we have to use expansive nonlinearities (exponential, powerlaw) sparingly
% because they can explode and get us into bad numerical territory

% specify the nonlinearities we'll fit with
f = @nlfuns.logexp1;
g = {@nlfuns.threshLinear, @nlfuns.logistic};

% regularizaiton (this imposes smoothness in 2D space)
Cinv = qfsmooth2nd(nx, ny); % 2nd derivative matrix
Cinv = Cinv + eye(nx*ny); % Identity matrix --> ridge regression

% start by fitting the LNP model
mstruct = struct();
mstruct.neglogli = @regression.neglogli_poiss; % set the likelihood here
mstruct.logprior = @regression.logprior_Cinv;  % set the prior here
Xd = [X ones(size(X,1),1)]; % augment stimulus with a column of ones -- this creates a bias term

inds = find(targetCond==1);
mstruct.liargs = {Xd(inds,:), R(inds), f, binSize};
mstruct.priargs = {blkdiag(Cinv, 0)};
fun = @(w) regression.neglogpost_GLM(w, 10, mstruct);

% options for fitting (turn on gradient checking)
% opts = optimset('GradObj', 'On', 'Hessian', 'Off', 'Display', 'none', 'MaxIter', 10e10);
% opts = optimoptions('fminunc','Algorithm','quasi-newton','Display', 'Final', 'SpecifyObjectiveGradient',true, 'MaxIter', 10e10);
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display', 'iter', 'SpecifyObjectiveGradient',true, 'MaxIter', 10e10, 'CheckGradients', true, 'OptimalityTolerance', 1e-9, 'FiniteDifferenceType', 'central');
w1 = fminunc(fun, [w0; 0], opts);

figure(1); clf
subplot(1,3,1)
imagesc(reshape(w0, [ny nx]))
title('Spike Triggred Average')

subplot(1,3,2)
imagesc(reshape(w1(1:end-1), [ny nx]))
title('LNP')

figure(2);
subplot(2,1,1); cla
plot(R, 'k');
hold on
RhatLN = f(Xd*w1)*binSize;
plot(RhatLN, 'r')
xlabel('Sample #')
ylabel('Spike Count')
legend('Data', 'LNP')

for i = nTargets
    inds = find(targetCond==i);
    subplot(2,nConds,nConds+find(i==nTargets))
    sta = simpleSTC(X(inds,:),RhatLN(inds), 1);
    sta = reshape(sta, [ny nx]);
    imagesc(sta)
    title(sprintf('STA (%d targs)', i))
    
end
set(gcf, 'PaperSize', [6 3], 'PaperPosition', [0 0 6 3])
saveas(gcf, fullfile(figDir, 'lnpSTA.png'))


% initialize the gain field by fitting a hacky modulated poisson
Xd = [RhatLN(:) bsxfun(@times, X, RhatLN(:)) ones(numel(R),1)];

mstruct.liargs = {Xd(inds,:), R(inds), f, binSize};
mstruct.priargs = {blkdiag(.1, Cinv, 0)};

fun = @(w) regression.neglogpost_GLM(w, 1, mstruct);
wmod = fminunc(fun, [1; w1], opts);

figure(1)
subplot(2,3,6)
imagesc(reshape(wmod(2:end-1), [ny nx]))
title('Gain field Init')

% initialize filters
wGainInit =  wmod(2:end-1);
wRFInit = w1;
%% Fit gain model
% change the likelihood function to the modulated poisson model
mstruct.neglogli = @regression.neglogli_modulated_poiss;

% likelihood arguments: 
% {stimulus, spikes, subunit nonlinearity, spike, nonlinearity, bin size}
Xd = [X ones(size(X,1),1)]; % back to the original stimulus
inds = find(targetCond==1); % fit the RF using the one-target case
mstruct.liargs = {Xd(inds,:), R(inds), g{1}, f, binSize};
mstruct.priargs = {blkdiag(Cinv, 0.1)};

maxIter = 100; % maximum iterations to try

% set up the cost function
fun = @(w) regression.neglogpost_GLM(w, rho, mstruct);

wInitRF = w1; 
% wInitRF = max(w1,0); % initialize with LNP, but only take positive values

% trying other initialization options
% wHatRF = wInitRF;
% wHatRF = [Ke(:); rfOffset];
wHatRF = fminunc(fun, wInitRF, opts);

figure(1)
subplot(2,3,3)
imagesc(reshape(wHatRF(1:end-1), [ny nx]))
title('RF')

% initialize denominator

% then generate the output of the RF, that gets passed into the likelihood
% function for the modulated poisson model
numerator = g{1}(Xd*wHatRF); % evaluate the signal in the numerator
% numerator = f(Xd*w1); % use the LNP model as the modulator signal

% swap out likelihood arguments for the divisive term 
mstruct.liargs = {Xd, R, g{2}, f, binSize, numerator};
fun = @(w) regression.neglogpost_GLM(w, rho_sup, mstruct);

% fun = @(w) regression.neglogli_modulated_poiss(w, X, R, g{2}, f, binSize, numerator);
% initialize with RF estimate? use STC?
% wInitGain = [Ks(:); dfOffset]; % if initializing with true filter
wInitGain = [wGainInit(:); 0.1]; %rand(size(wHatRF));
% wInitGain = rand(numel(wHatRF),1); % random initialization

[wHatGain, fInit] = fminunc(fun, wInitGain, opts);

figure(1)
subplot(2,3,3)
imagesc(reshape(wHatRF(1:end-1), [ny nx]))
title('RF')

subplot(2,3,6)
imagesc(reshape(wHatGain(1:end-1), [ny nx]))
title('div RF')
drawnow

% turn off gradient checking -- it fails sometimes  (I believe there are
% numerical errors from numerical chain rule the way I implemented it). If
% this continues
opts = optimoptions('fminunc','Algorithm','quasi-newton','Display', 'None', 'SpecifyObjectiveGradient',true, 'MaxIter', 10e10, 'OptimalityTolerance', 1e-9);

modulator = g{2}(Xd*wHatGain);

figure(2); 
subplot(2,1,1)
plot(f(numerator.*modulator)*binSize, 'g')

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
    mstruct.priargs = {blkdiag(Cinv, 1e3)};
    fun = @(w) regression.neglogpost_GLM(w, rho, mstruct);
    
    [wHatRF,~] = fminunc(fun, wInitRF, opts);
    
    subplot(1,2,1)
    plot(wHatRF);
    title('RF')
    
    % fit denominator
    numerator = g{1}(Xd*wHatRF);
    mstruct.liargs = {Xd, R, g{2}, f, binSize, numerator};
    fun = @(w) regression.neglogpost_GLM(w, rho_sup, mstruct);

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
subplot(2,2,1)
imagesc(Ke)
title('K_e True')
subplot(2,2,3)
imagesc(reshape(wHatRF(1:end-1), [ny nx]))
title('K_e fit')

subplot(2,2,2)
imagesc(Ks)
title('K_s True')

subplot(2,2,4)
imagesc(reshape(wHatGain(1:end-1), [ny nx]))
title('K_s fit')

set(gcf, 'PaperSize', [5 5], 'PaperPosition', [0 0 5 5])
saveas(gcf, fullfile(figDir, 'divSfits.png'))

Rhat = f(g{1}(Xd*wHatRF).*g{2}(Xd*wHatGain))*binSize;

figure(2); clf
subplot(2,1,1)
plot(R, 'k'); hold on
plot(Rhat, 'g')
xlabel('Sample #')
ylabel('Spike Count')
legend('Data', 'divS')

for i = nTargets
    inds = find(targetCond==i);
    subplot(2,nConds,nConds+find(i==nTargets))
    sta = (X(inds,:)'*R(inds))/sum(R(inds));
    sta = reshape(sta, [ny nx]);
    imagesc(sta, [0 .09]);
    
end


set(gcf, 'PaperSize', [6 3], 'PaperPosition', [0 0 6 3])
saveas(gcf, fullfile(figDir, 'divSSTA.png'))


