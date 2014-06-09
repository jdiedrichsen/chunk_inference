function [chunk_means_rt, rt_cov, chunk_means_er, er_cov] = ...
    create_chunk_means_covs(chunks, cor_chunks, ...
        mean_pause, mean_inchunk, v, rho, ...
        mean_pause_er, mean_inchunk_er, v_er, rho_er)

n_chunks = size(chunks, 1);

% Correlation
rt_cov = zeros(size(cor_chunks));
er_cov = zeros(size(cor_chunks));
n_seq_len = size(cor_chunks, 1);
for z = 1:n_chunks
    rt_cov(:, :, z) = ((cor_chunks(:, :, z) - eye(n_seq_len))*rho + ...
        eye(n_seq_len))*v;
    er_cov(:, :, z) = ((cor_chunks(:, :, z) - eye(n_seq_len))*rho_er + ...
        eye(n_seq_len))*v_er;
end

% mean of chunks + mean at beginning of chunk
chunk_means_rt = ...
    (diff([zeros(size(chunks, 1), 1) chunks], 1, 2)>0)*mean_pause + ...
    (~diff([zeros(size(chunks, 1), 1) chunks], 1, 2)>0)*mean_inchunk;
chunk_means_er = ...
    (diff([zeros(size(chunks, 1), 1) chunks], 1, 2)>0)*mean_pause_er + ...
    (~diff([zeros(size(chunks, 1), 1) chunks], 1, 2)>0)*mean_inchunk_er;