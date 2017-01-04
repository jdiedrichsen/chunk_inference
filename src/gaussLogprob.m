function logp = gaussLogprob(arg1, arg2, arg3)
% Log pdf of multivariate Gaussian
% L  = gaussLogprob(model, X)
% L = gaussLogprob(mu, Sigma, X)
% X(i,:) is i'th case, can contain NaNs for missing values.
% Sigma can be a vector; this is interpreted as a diagonal matrix.
% *** In the univariate case, Sigma is the variance, not the standard
% deviation! ***
%
% Examples
% L  = gaussLogprob(zeros(1,3), randpd(3), rand(10,3))
% L  = gaussLogprob(zeros(3,1), diag(randpd(3)), rand(10,3))
% L = gaussLogprob(structure(mu, Sigma), X)
%
% Change from original: The mean vector is supposed to be same orientation
% as the data - that is it's a row vector for constant mean and a matrix of
% the ssame size as the data for variable mean
% This file is from pmtk3.googlecode.com


switch nargin
    case 3,  mu = arg1; Sigma = arg2; X = arg3;
    case 2, model = arg1; mu = model.mu; Sigma = model.Sigma; X = arg2;
    otherwise
        error('bad num args')
end

% Deal with different size of the mu argument
if isscalar(mu)
    X = X(:);
end
[N, d] = size(X);
[d1,d2] = size(mu);
if (d2~=d) 
    error('mu must be either scalar of have the same number of columns as the data'); 
end; 

switch(d1) 
    case 1
        X = bsxfun(@minus, X, rowvec(mu));
    case N 
        X = X-mu; 
    otherwise 
        error('mu needs to have either 1 row (constant mean) or match the number of rows in the data'); 
end; 

if any(isnan(X(:)))
    logp = gaussLogprobMissingData(mu,Sigma, X);
    return;
end

if isvector(Sigma) && (numel(Sigma) > 1) % diagonal case
    sig2 = repmat(Sigma', N, 1);
    tmp  = -(X.^2)./(2*sig2) - 0.5*log(2*pi*sig2);
    logp = sum(tmp, 2);
else
    % Full covariance case
    if 0
        logp2 = -0.5*sum((X/Sigma).*X, 2);
        logZ2 = (d/2)*log(2*pi) + 0.5*logdet(Sigma);
        logp2 = logp2 - logZ2;
    end
    %   slightly faster version
    R    = chol(Sigma);
    logp = -0.5*sum((X/R).^2, 2);
    logZ = 0.5*d*log(2*pi) + sum(log(diag(R)));
    logp = logp - logZ;
    %assert(approxeq(logp, logp2))
end

end

