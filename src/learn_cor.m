function rho = learn_cor(chunks, cor_chunks, data, gamma, ...
    start_pause, nonstart_pause, variance, fit_mean)
% learn correlation

ind_chunk_start = diff([zeros(size(chunks, 1), 1) ...
    chunks], 1, 2)>0;

n_numerator = 0;
n_chunks = size(chunks, 1);
n_seq_len = size(chunks, 2);
tmp_cov = 0;
tmp_var = 0;
n_COV   = 0; 
n_VAR   = 0; 
goodData=~isnan(data); 
for i = 1:n_chunks
    if fit_mean
        % center movement_times
        cmt = bsxfun(@minus, data, ...
            ind_chunk_start(i, :)*start_pause + ...
            (~ind_chunk_start(i, :))*nonstart_pause);
    else
        cmt = mt_seq;
    end
    % Old version: 
    % Sigma = ((bsxfun(@times, cmt, gamma(2:end, i)))'*cmt) .* ...
    %     cor_chunks(:, :, i) .* (1-eye(n_seq_len));
    % n_numerator = n_numerator + ...
    %     sum(gamma(2:end, i))*nnz(cor_chunks(:, :, i) .* (1-eye(n_seq_len)));

    % To ignore the NaN observations: 
    
    % Places in the covariance matrix we want to look at: 
    indCOV = cor_chunks(:,:,i) .* (1-eye(n_seq_len)); 
    % Corresponding places in the variance matrix we want to look at 
    indVAR = diag(sum(indCOV,2)); 
    
    
    SS  = bsxfun(@times,permute(cmt,[1 2 3]),permute(cmt,[1 3 2]));  % Sums of squares for each individual trial
    wSS = bsxfun(@times,SS,gamma(:, i));    % Weighted by the posterior probaility 
    sSS = squeeze(nansum(wSS,1));               % Take sum ignoring the NaNs
    COV = sSS .* indCOV;
    VAR = sSS .* indVAR; 
    num = (bsxfun(@times, goodData, gamma(:, i)))'*goodData;
    numCOV   = num .* indCOV;  
    numVAR   = num .* indVAR; 
    tmp_cov     = tmp_cov + nansum(COV(:));
    tmp_var     = tmp_var + nansum(VAR(:)); 
    n_COV    = n_COV + sum(numCOV(:)); 
    n_VAR    = n_VAR + sum(numVAR(:)); 
end
var2    = tmp_var./n_VAR;
new_cov = tmp_cov./n_COV;  
rho = new_cov/var2;
if (abs(rho)>1)
    keyboard; 
end; 