clear all, clc, close all
%% Load/Generate Data

% Y is an n x p <n samples  x  p features> matrix
% load(...)

% DEMO: generate data (spiked cov model)
p = 1000;                   % ambient dimension
n = 400;                    % number of samples
k = 50;                     % cardinality of support
support = 1:k;              % an arbitrary support of cardinality k
snr = 20;

x = sparse(support, 1, ...  % random sparse signal
           rand(k, 1),...    
                p, 1);      
x = x / norm(x);

Y = sqrt(snr) * diag(randn(n, 1)) * repmat(x', n, 1) + randn(n, p); % Data


%% Sparsity Constraint

% Set up
params.algorithm     = 'sparse';          % Seek sparse solution
params.nnz           = [10, 30,  50, 70]; % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 1e4;               % Number of subspace samples

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = '';                % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window


% Run
[X] = spanpc(Y, params);

% Plot explained variance for each sparsity value
figure;
plot(params.nnz, var(Y*X), '--sr');
title('Explained variance: k-sparse principal component')
xlabel('Sparsity (k)');
ylabel('Explained (empirical) Variance');
grid on;

clear params

%% Nonnegativity & Sparsity Constraints

% Set up
params.algorithm     = 'nnsparse';        % Seek sparse solution
params.nnz           = [10, 30,  50, 70]; % Target sparsity values
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 1e4;               % Number of subspace samples

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = 'sample.log';      % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window


% Run
[X] = spanpc(Y, params);

% Plot explained variance for each sparsity value
figure;
plot(params.nnz, var(Y*X),'--sr');
title('Explained variance: nonnegative k-sparse principal component')
xlabel('Sparsity (k)');
ylabel('Explained (empirical) Variance');
grid on;

clear params

%% Nonnegativity Constraint (No sparsity)

% Set up
params.algorithm     = 'nonnegative';     % Seek sparse solution
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 1e4;               % Number of subspace samples

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = '';                % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'off';             % Log level for command window

% Run
[x] = spanpc(Y, params);

% Explained variance for each sparsity value
var(Y*x)

clear params

%% ST Path on a Graph

% adj = [
%      0     0     0     0     1     1     1     1     0     0     0     0
%      0     0     0     0     1     1     1     1     0     0     0     0
%      0     0     0     0     1     1     1     1     0     0     0     0
%      0     0     0     0     1     1     1     1     0     0     0     0
%      0     0     0     0     0     0     0     0     1     1     1     1
%      0     0     0     0     0     0     0     0     1     1     1     1
%      0     0     0     0     0     0     0     0     1     1     1     1
%      0     0     0     0     0     0     0     0     1     1     1     1
%      0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0
%      0     0     0     0     0     0     0     0     0     0     0     0
%      ];
% S = [1,  2,  3,  4]; %<- all these vertices as sources
% T = [9, 10, 11, 12]; %<- all these vertices as terminals
% Gen graph
graphparam.n = p;
graphparam.L = round(p/100);
[adj, layerSizes] = generateTrellisGraph(graphparam);
S = 1:layerSizes(1);
T = p-layerSizes(end)+1:p;

% -----------------------------------------------
% Here you put your input
Gst = graphaddst(adj, S, T);
spy(Gst.adj)
% Set up
params.algorithm     = 'stpath';          % Seek sparse solution
params.apprxrank     = 3;                 % Rank of approximation
params.maxsamples    = 1e4;               % Number of subspace samples
params.graph         = Gst;

% (Optional)
params.inputdata     = 'rows';            % Samples in Y are rows
params.maxnoupditer  = 1e3;               % Max iters without improvement    
params.centerdata    = true;              % Center samples (subtract mean)
params.standardata   = false;             % Standardize features
params.logfile       = '';                % Path to log file
params.logfilelevel  = 'off';             % Log level for log file
params.logcwlevel    = 'info';            % Log level for command window

% Run
[x] = spanpc(Y, params);

% Explained variance for each sparsity value
var(Y*x)

clear params