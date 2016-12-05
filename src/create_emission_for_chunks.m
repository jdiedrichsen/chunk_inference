function p_obs_un = create_emission_for_chunks(P,Data,Model); 
% Create emission structure for chunks.
% Create probability of observation and movement time covariances 
% according to chunk.
%
% Parameters:
%    rho        : correlation within chunks
%    rho_er     : correlation of errors
%    self_t     : self transition probability
%    mt_seq     : detrended movement times observed
%    er_seq     : detrended error sequence
%    mean_pause : pause at beginning of chunk
%    
% Returns:
%    mt_cov: movement time covariance for all chunks
%    p_obs : probability of observation
%    T     : transition probabilities across chunks


if ~Model.fit_er && ~Model.fit_rt
    error('Either fit_er or fit_mt must be true');
end
% Find best parameters for chunking

if isempty(Model.chunks)
    chunks = create_chunks3();    
end
cor_chunks = to_corr_chunks(Model.chunks);
n_chunks = size(Model.chunks, 1);

if (Model.fit_rt) 
    [chunk_means_mt, mt_cov] = create_chunk_means_covs(Model,P,'mt'); 
end; 
if (Model.fit_er) 
    [chunk_means_er, er_cov] = create_chunk_means_covs(Model,P,'er'); 
end; 

% Compute probability of observation
n_time = size(Data.mt_seq, 1);
p_obs = zeros(n_time, n_chunks);

for i = 1:n_chunks
    if Model.fit_rt && Model.fit_er    
        p_obs(:, i) = ...
            gaussLogprob(chunk_means_mt(i, :)', ...
                    mt_cov(:, :, i), Data.mt_seq) + ...
            gaussLogprob(chunk_means_er(i, :)', ...
                    er_cov(:, :, i), Data.er_seq);
    elseif Model.fit_rt
        p_obs(:, i) = ...
            gaussLogprob(chunk_means_mt(i, :)', ...
            mt_cov(:, :, i), Data.mt_seq);
    else
        p_obs(:, i) = gaussLogprob(chunk_means_er(i, :)', ...
            er_cov(:, :, i), Data.er_seq);
    end
end
% p_obs = exp(normalizeLogspace(p_obs));
p_obs = exp(p_obs);
% Unnormalized
p_obs_un = p_obs;

% Normalize
% p_obs = bsxfun(@rdivide, p_obs, sum(p_obs, 2));
