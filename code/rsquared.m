function r2 = rsquared(rtrue, rhat, rbase)
% r2 = rsquared(rtrue, rhat)
if nargin < 3
    rbase = mean(rtrue(:));
end

r2 = 1-(sum((rtrue(:)-rhat(:)).^2))/(sum((rtrue(:)-rbase).^2));