function ppca(T)
% A script that performs PPCA


% E-step

% M-step


% Auxiliary functions

% Bound
function bound(T)
Lc = 0;
[N,D] = size(T);

for it = 1:N
  Lc = 0.5*d*log(sig2) + 0.5*trace(exxt)+0.5*(sig2)^(-1)*(t(it,:)-mu(it,:))'*(t(it,:)-mu(it,:)) \
	  -(sig2)^(-1)*exn'*W'*(t(it,:)-mu(it,:))+0.5*(sig2)^(-1)*trace(W'*W*exxt)
end

% Posterior moments
function [ex, exxt] = computeMoments(t)

invM = inv(M);
exn = invM*W'*(t-mu);
exxt = sig2*invM + exxt;