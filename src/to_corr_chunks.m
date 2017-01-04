function cor_chunks = to_corr_chunks(chunks,cov_structure)
% Transform chunks into predicted correlation structure 
% cor_chunks = to_corr_chunks(chunks,cov_structure)
cor_chunks = zeros(size(chunks, 2), size(chunks, 2), size(chunks, 1));
switch cov_structure
    case 'none' 
         for z = 1:size(chunks, 1)
                cor_chunks(:,:,z) = eye(size(chunks,2)); 
         end; 
    case 'chunk'
        for z = 1:size(chunks, 1)
            cur_d = 0;
            for j = 1:size(chunks, 2)
                if chunks(z, j) == 0
                    cur_d = 0;
                    cor_chunks(j, j, z) = 1;
                elseif (j>1 && chunks(z, j) ~= chunks(z, j-1))
                    cur_d = 1;
                    cor_chunks(j, j, z) = 1;
                else
                    cor_chunks(max(j-cur_d, 1):j, ...
                        max(j-cur_d, 1):j, z) = 1;
                    cur_d = cur_d + 1;
                end
            end
        end
        
    case 'category'
        cor_chunks = NaN(size(chunks, 2), size(chunks, 2), size(chunks, 1));
        
        diffchunks = diff(chunks,1,2);
        diffchunksnormal = [ones(size(diffchunks,1),1) diffchunks];
        diffchunksinv = [zeros(size(diffchunks,1),1) ~diffchunks];
        for p = 1:size(chunks,1)
            x = find(diffchunksinv (p,:)== 0);
            for i = 1:size(chunks,2)
                if any(ismember(x,i))
                    cor_chunks(i,:,p) = diffchunksnormal(p,:);
                else
                    cor_chunks(i,:,p) = diffchunksinv(p,:);
                    
                end
            end
        end
        
end