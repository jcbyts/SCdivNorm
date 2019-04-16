%%
fun = @(w) regression.neglogli_modulated_poiss(w, X, R, g{2}, f, .1, numerator);

w0 = rand(size(X,2),1);
[f0, df] = fun(w0);

rvec = randn(size(X,2),1)*10e-6;

df0 = df'*rvec;
f1 = fun(w0+rvec);

[f1-f0 df0]
%%

% fun = @(w) regression.neglogli_poiss(w, X, R, f, binSize);
for i = 1:1000
[d1(i), d2(i)] = DerivCheck(fun, rand(size(X,2),1), opts);
end

figure(2); clf
plot(d1, d2, '.')
hold on
plot(xlim, xlim, 'k')

figure(3); clf
plot(d1-d2, '.')