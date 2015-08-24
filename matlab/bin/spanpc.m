function [X] = spanpc(data, params)

params= spanpc_parse_input(data, params);

if strcmpi(params.inputdata, 'columns')
    data = data';
end

% s : number of samples (rows of X)
% p : number of features (columns of X)
[s, p] = size(data);

LOG = params.xx_logger;


LOG.debug(mfilename, evalc('disp(params)'));

% Preprocess data
if params.centerdata              % Remove mean from samples
    samplemean = mean(data, 1);
    data = data - repmat(samplemean, s, 1);
end
if params.standardata              % Normalize nonzero features
    data = data * diag(mean(data.^2, 1).^0.5)^-1;
end


% Adjust sparsity params if such exists
if isfield(params, 'xx_nnz')
   params.xx_nnz = min(params.xx_nnz, p);
end

LOG.info(mfilename, sprintf('Computing principal subspace...'));
r = params.apprxrank;

% Compute V <p x r>: r leading components of S
[~, Sx, Vx] = svds(data, r);
V = sqrt(s)^(-1) * Vx * Sx;  % COV(X) = (1/S) X'*X = Vx*Sx*Sx'*Vx';

X = subspacesampler(V, params);

end % end of function

%--------------------------------------------------------------------------
function [X] = subspacesampler(V, params)

LOG         = params.xx_logger;
rank1solver = params.xx_rank1solver;

% Constants
SMPL_TOTAL  = params.maxsamples;
BATCH_SIZE  = params.xx_samplesperbatch;
SMPL_EARLY  = params.xx_earlystopmax;

% Initialization
[p, d]      = size(V);
smpl_rem    = SMPL_TOTAL;
smpl_batch  = params.xx_samplesfirstbatch;
num_cand    = params.xx_candpersample;
scores      = zeros(1, num_cand);
X           = sparse([], [], [], p, num_cand, 0);

flg_earlystop = false;
last_upd    = 0;

% Run
LOG.info(mfilename, sprintf('Subspace sampling...'));
T = tic;

while(true)
    
    batch_start = tic;
    
    C = randn(d, smpl_batch);
    for i = 1:smpl_batch

        c = C(:,i) / norm(C(:,i));
        w = V*c;
        X_ = rank1solver(w, params);

        for cand = 1:num_cand 
            
            score_ = norm(V'*X_(:, cand));
            if score_ > scores(cand)
                X(:, cand)   = X_(:, cand);
                scores(cand) = score_;
                last_upd = 0;
            end
            
        end % updated all candidates
        
        last_upd = last_upd + 1;
        if last_upd >= SMPL_EARLY
            flg_earlystop = true;
            break;
        end
        
    end
    
    batch_time = toc(batch_start);
    
    if flg_earlystop
        LOG.info(mfilename, ... 
            sprintf('No update in %d iter (Early stop)', SMPL_EARLY));
        break; 
    end
    
    smpl_rem = smpl_rem - smpl_batch;
    if smpl_rem == 0
        break;
    end
    
    LOG.info(mfilename, ...
             sprintf('%3.0f%% in %d seconds (%d left)', ...
                      (1 - smpl_rem/SMPL_TOTAL)*100, ...
                      round(toc(T)), ...
                      round(smpl_rem * batch_time / smpl_batch)) );
    
    smpl_batch = min([BATCH_SIZE, smpl_rem]);
end
LOG.info(mfilename, 'Completed.');


end % end of function

