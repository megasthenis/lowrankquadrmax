function [X] = rank1solver_nnsparse(w, options)
% RANK1SOLVER_NNSPARSE Solve rank-1 nonnnegative sparse PCA problem:
%                    min (w'x)^2, 
%             subject to norm(x,2) = 1, norm(x,0) <= k, x >= 0
% 
% [X] = rank1solver_sparse(w, options) compute the solution for argument w
% and sparsity k determined in the options.nnz parameter. 
% X has number of columns equal to numel(options.nnz). Each column of X 
% contains the solution corresponding to a different sparsity value.

p = numel(w);

pos = sum(w >= 0);                         % num of positive entries
neg = p - pos;                             % num of negative entries
    
[~, indxsort] = sort(w, 'descend');

%A candidate for each sparsity value
X = sparse([], [] ,[], p, options.xx_candpersample, 0);

for cand = 1:options.xx_candpersample

    k = options.xx_nnz(cand);              % target sparsity

    k_pos = min([k, pos]);                 % effective sparsity (pos case)
    spprt_pos = indxsort(1:k_pos);         % candidate support  (pos case)
    val_pos = norm(w(spprt_pos), 2);

    k_neg = min([k, neg]);                 % effective sparsity (neg case)
    spprt_neg = indxsort(end+1-k_neg:end); % candidate support  (neg case)
    val_neg = norm(w(spprt_neg), 2);

    if val_pos >= val_neg
        X(spprt_pos, cand) = + w(spprt_pos) / val_pos;
    else
        X(spprt_neg, cand) = - w(spprt_neg) / val_neg;

    end

end

end % end of function