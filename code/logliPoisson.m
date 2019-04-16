function llsp = logliPoisson(lambda, r)
% Poisson log likelihood

etol = 1e-100;
lambda(lambda<etol)=etol;

llsp   = r'*log(lambda) - sum(lambda);