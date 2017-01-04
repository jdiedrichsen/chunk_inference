function [M, COV] = create_chunk_means_covs(Model, P, type, trialNum)

n_chunks = size(Model.chunks, 1);
n_seq_len = size(Model.cor_chunks, 1);

switch(type)
    case 'mt'
        boundary = diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0;
        COV = zeros(size(Model.cor_chunks));
        for z = 1:n_chunks
            COV(:, :, z) = ((Model.cor_chunks(:, :, z) - eye(n_seq_len))*P.rho + eye(n_seq_len));
            if (Model.fit_rt_var==1)
                COV(:,:,z) = COV(:,:,z)*P.v;
            elseif(Model.fit_rt_var==2)
                v=sqrt(boundary(z,:)*P.v(1)+(1-boundary(z,:))*P.v(2));
                COV(:,:,z)=COV(:,:,z).*(v'*v);
            end;
        end;
        boundary = permute(boundary,[3 2 1]);
        if (Model.fit_rt==1)
            M = boundary*P.mean_pause + (~boundary)*P.mean_inchunk;
        elseif (Model.fit_rt==2)
            pred_pause  =  (P.pause_B-P.pause_A)*exp(-trialNum/P.pause_tau)+P.pause_A;
            pred_inchunk = (P.inchunk_B-P.inchunk_A)*exp(-trialNum/P.inchunk_tau)+P.inchunk_A;
            M = bsxfun(@times,pred_pause,boundary)+bsxfun(@times,pred_inchunk,~boundary);
        end;
    case 'er'
        COV = zeros(size(Model.cor_chunks));
        for z = 1:n_chunks
            COV(:, :, z) = ((Model.cor_chunks(:, :, z) - eye(n_seq_len))*P.rho_er + ...
                eye(n_seq_len))*P.v_er;
        end;
        boundary = diff([zeros(size(Model.chunks, 1), 1) Model.chunks], 1, 2)>0;
        boundary = permute(boundary,[3 2 1]);
        M = boundary*P.mean_pause + (~boundary)*P.mean_inchunk;
end;



