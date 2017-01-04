function logp = gaussLogprobMissingData(arg1, arg2, arg3)
% Same as gaussLogprob, but supports missing data, represented as NaN
% X is N*D
% L = gaussLogprobMissingData(mu, Sigma, X)
% L = gaussLogprobMissingData(model, X)

% Test example
%{
N = 5; D = 3;
X = randn(N,D);
M = rand(N,D)>0.5;
X(M) = NaN;
mu = zeros(D,1); Sigma = eye(D);
model = gaussCreate(mu, Sigma);
logp = gaussLogprobMissingData(model, X);
logp2 = gaussLogprobMissingData(mu, Sigma, X);
assert(approxeq(logp, logp2))
%}

% This file is from pmtk3.googlecode.com

switch nargin
    case 3,  mu = arg1; Sigma = arg2; X = arg3;
    case 2, model = arg1; mu = model.mu; Sigma = model.Sigma; X = arg2;
    otherwise
        error('bad num args')
end

[n,d] = size(X); 




logp = NaN(n,1);
missRows = ~isnan(X);

% Find the patterns of missing data -
% Then apply the multivaraite Gaussian to each pattern seperately
[missPattern,~,pat]=unique(missRows,'rows');
numPatterns = size(missPattern,1);
for i=1:numPatterns
    indx = find(pat==i);
    nnan = missPattern(i,:);
    logp(indx) = gLP(Sigma(nnan,nnan), X(indx,nnan));
end;

function logp = gLP(Sigma,X);
d = size(X,2); 
R    = chol(Sigma);
logp = -0.5*sum((X/R).^2, 2);
logZ = 0.5*d*log(2*pi) + sum(log(diag(R)));
logp = logp - logZ;
