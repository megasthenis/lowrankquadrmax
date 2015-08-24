function [X] = rank1solver_sparse(w, options)
% RANK1SOLVER_SPARSE Solve rank-1 sparse PCA problem:
%                    min (w'x)^2
%             subject to norm(x,2) = 1, norm(x,0) <= k
% 
% [X] = rank1solver_sparse(w, options) compute the solution for argument w
% and sparsity k determined in the options.nnz parameter. 
% X has number of columns equal to numel(options.nnz). Each column of X 
% contains the solution corresponding to a different sparsity value.

p = numel(w);

[~, ind_order] = sort(abs(w), 'descend');

%A candidate for each sparsity value
X = sparse([], [] ,[], p, options.xx_candpersample, 0);

for cand = 1:options.xx_candpersample
    
    k = options.xx_nnz(cand);                 % target sparsity

    support = ind_order(1:k);

    X(support, cand) = w(support) / norm(w(support), 2);

end

end % end of function