function [M, COV] = create_chunk_means_covs(Model, P, type)

n_chunks = size(Model.chunks, 1);
n_seq_len = size(Model.cor_chunks, 1);

switch(type)
    case 'mt'
        COV = zeros(size(Model.cor_chunks));
        for z = 1:n_chunks
            COV(:, :, z) = ((Model.cor_chunks(:, :, z) - eye(n_seq_len))*P.rho + ...
                eye(n_seq_len))*P.v;
        end;
        M = (diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0)*P.mean_pause + ...
            (~diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0)*P.mean_inchunk;
    case 'er'
        COV = zeros(size(Model.cor_chunks));
        for z = 1:n_chunks
            COV(:, :, z) = ((Model.cor_chunks(:, :, z) - eye(n_seq_len))*P.rho_er + ...
                eye(n_seq_len))*P.v_er;
        end;
        M = (diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0)*P.mean_pause_er + ...
            (~diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0)*P.mean_inchunk_er;
end;



