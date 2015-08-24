function [ params ] = spanpc_parse_input(data, params)

[s, p] = size(data);

PARAMLIST = {...
    %{ 
        parameter name,
        default value, 
        required?,  (if false then should have default value)
        valid values, 
        description,
        validating function
    %}
    {
        'algorithm',...
        [],...
        true,... 
        {'(String) sparse|nnsparse|nonnegative'},...
        'Specifies the type of constraints.',...
        @(val) instringset(val, {'sparse', 'nnsparse', 'nonnegative', 'stpath'}), ...
    };
    {
        'inputdata',...
        'rows',...
        false,... 
        {'(String) rows (default)|columns|'},...
        ['Specifies if the samples are the rows or columns of the ',...
         'input data matrix.'
        ],...
        @(val) instringset(val, {'rows', 'columns'}), ...
    };
    {
        'apprxrank', ...
        [],...
        true,...
        {'(Positive Integer)'}, ...
        ['The rank of the low-rank approximation to be used.',... 
         'Higher values lead to better results, but longer execution time.'
        ], ...
        @(val) isposint(val), ...
    };
    {
        'nnz',...
        [],...
        false,...
        {'(Positive Integer)'}, ...
        ['Controls the sparsity (number of non-zero entries) of the ', ...
         'output. If nnz contains multiple values, ', ...
         'an output is computed for each value.'
        ],...
        @(val) isposintorinfarray(val), ...
    };
    {
        'centerdata',...
        false,...
        false,...
        {'true|false'}, ...
        ['Compute and remove the sample mean from the data',...
         ' prior to applying the algorithm.'
        ], ...
        @(val) islogical(val), ...
    };
    {
        'standardata',...
        false,...
        false,...
        {'true|false'}, ...
        ['Normalize each feature with its sample variance',...
         ' prior to applying the algorithm on the data.'
        ], ...
        @(val) islogical(val), ...
    };
    { 
        'maxsamples', ...
        [],...
        true,...
        {'(Positive Integer)'},...
        ['The number of samples to consider from the principal subpsace.' ...
         'More samples imply possibly improved accuary, ',...
         'but longer execution time.'
        ],...
        @(val) isposint(val), ...
    };
    { 
        'logfilelevel', ...
        'off',...
        false,...
        {'(String) all|trace|debug|info|warn|error|fatal|off'},...
        ['Level of log detail for log file.'],...
        @(val) instringset(val, {'all', 'trace', 'debug', 'info',...
                                 'warn', 'error','fatal', 'off', ''} ...
                ), ...
    };
    { 
        'logfile', ...
        [],...
        false,...
        {'(String filename)'},...
        ['Path to log file.'],...
        @(val) isempty(val) || ischar(val), ...
    };
    { 
        'logcwlevel', ...
        'off',...
        false,...
        {'(String) all|trace|debug|info|warn|error|fatal|off'},...
        ['Level of log detail for command window.'],...
        @(val) instringset(val, {'all', 'trace', 'debug', 'info',...
                                 'warn', 'error','fatal', 'off', ''} ...
                ), ...
    };
    {
        'maxnoupditer', ...
        Inf,...
        false,...
        {'(Positive Integer)'},...
        ['Max number of samples to be considered without improving', ...
         'the objective function. If no improvement is observed ',...
         'after maxnoupditer samples, execution terminates.'
        ],...
        @(val) isinf(val) || isposint(val), ...
    };
};

            
% General check
for i = 1:size(PARAMLIST, 1)
    
    param_name      = PARAMLIST{i}{1};
    param_default   = PARAMLIST{i}{2};
    param_req       = PARAMLIST{i}{3};
    param_validator = PARAMLIST{i}{6};
   
    if ~isfield(params, param_name)
        if param_req
            error('required option %s not specified.', param_name);
        end
        
        % Parameter not required and not specified. Apply default value.
        params.(param_name) = param_default;
    end
    
    param_value = params.(param_name);
    if ~param_validator(param_value)
        error(['invalid value (', param_value, ')',...
               ' for option ', param_name,...
              ] ...
        );
    end
    
end

% Check parameters conditioned on algorithm
switch lower(params.algorithm)
    
    case 'sparse'
        
        if ~isfield(params, 'nnz') ...
            || ~isposintorinfarray(params.nnz) ...
            || isempty(params.nnz)
            error('algorithm [%s] requires option nnz', params.algorithm);
        end
        params.xx_rank1solver = @rank1solver_sparse;
        params.xx_nnz = params.nnz;
        params.xx_candpersample = numel(params.nnz);
        
    case 'nnsparse'
        
        if ~isfield(params, 'nnz') ...
            || ~isposintorinfarray(params.nnz) ...
            || isempty(params.nnz)
            error('algorithm [%s] requires option nnz', params.algorithm);
        end
        params.xx_rank1solver = @rank1solver_nnsparse;
        params.xx_nnz = params.nnz;
        params.xx_candpersample = numel(params.nnz);
        
    case 'nonnegative'
        
        if isfield(params, 'nnz') 
            if ~isempty(params.nnz)
               warning(['Ignoring params nnz.',...
                        ' Use algorithm nnsparse to enforce sparsity.',...
                       ]);
            end
        end
        params.xx_rank1solver = @rank1solver_nnsparse;
        params.xx_nnz = Inf;
        params.xx_candpersample = 1;
        
    case 'stpath'
        
        if ~isfield(params, 'graph') ...
            || isempty(params.graph)
            error('algorithm [%s] requires option graph', params.algorithm);
        end
        if  ~isvalidgraphwithst(params.graph)
            error('Invalid value for parameter: graph');
        end
        if params.graph.n-2 ~= p
           error('mismatch between data dimension and graph size') 
        end
        
        params.xx_rank1solver = @rank1solver_stpath;
        params.xx_nnz = Inf;
        params.xx_candpersample = 1;
        
        
    otherwise
        
        error('invalid algorithm: ''%s''', params.algorithm);
        
end
    
% Early Stop --------------------------------------------------------------
if isfield(params, 'maxnoupditer') && ~isempty(params.maxnoupditer)
    params.xx_earlystopmax = params.maxnoupditer;
else
    params.xx_earlystopmax = Inf;
end

% Samples per sample batch
params.xx_samplesfirstbatch = min([100,  params.maxsamples]);
params.xx_samplesperbatch   = ceil(params.maxsamples / 10);

% Early Stop (END) --------------------------------------------------------

% Logging -----------------------------------------------------------------
LOG = log4m.getLogger();

% Set up log file
level = getLogLevelCode(params.logfilelevel);
if isfield(params, 'logfile') ...
    && ~isempty(params.logfile) && ~strcmp(params.logfile, '')

    LOG.setFilename(params.logfile);
    LOG.setLogLevel(level);

else
    % No log file specified. If log level is set, ignore
    if level >= 0 && level ~= log4m.OFF
        warning(['Ignoring option logfilelevel.',...
                 ' Parameter logfile not specified.',...
                ]);
    end
    LOG.setLogLevel(log4m.OFF);
end

% Log level for command window
if isfield(params, 'logcwlevel')
    level = getLogLevelCode(params.logcwlevel);
    LOG.setCommandWindowLevel(level);
end

params.xx_logger = LOG;
% Logging (END) -----------------------------------------------------------


end % function



% Auxiliary functions
function code = getLogLevelCode(level)
    switch lower(level)
        case 'all'
            code = log4m.ALL;
        case 'trace'
            code = log4m.TRACE;            
        case 'debug'
            code = log4m.DEBUG;
        case 'info'
            code = log4m.INFO;
        case 'warn'
            code = log4m.WARN;
        case 'error'
            code = log4m.ERROR;
        case 'fatal'
            code = log4m.FATAL;
        case 'off'
            code = log4m.OFF;
        case '' 
            code = log4m.OFF;
        otherwise
            code = -1;
    end
end

function is = isposint(x)
    is =   numel(x) == 1 ...
        && x >= 0 ...
        && mod(x, 1) == 0;
end

function is = isposintarray(x)
    is = true;
    for i = 1:numel(x)
        if ~isposint(x(i))
           is = false;
           break;
        end
    end    
end

function is = isposintorinfarray(x)
    is = true;
    for i = 1:numel(x)
        if ~( isposint(x(i)) || isinf(x(i)) )
           is = false;
           break;
        end
    end    
end

function is = instringset(val, set)
    is = false;
    for i = 1:numel(set)
        if strcmpi(val, set{i})
           is = true;
           break;
        end
    end
end

function is = isvalidgraphwithst(G)
    is = 1;
    if ~isfield(G, 'adj') || isempty(G.adj) ...
        || any(G.adj(:) < 0) || any(G.adj(:) > 1) ...
        || any(mod(G.adj(:), 1) ~= 0)
        is = 0;
        return;
    end
    
    [n, junk] = size(G.adj);
    if n~= junk
       is = 0;
       return;
    end
    
    if ( ~isfield(G, 's') || isempty(G.s) ...
            || ~isposint(G.s) || (G.s > n) || (G.s < 1) )
        is = 0;
        return;
    end
    
    if ~isfield(G, 't') || isempty(G.t) ...
            || ~isposint(G.t) || G.t > n || G.t < 1
       is = 0;
       return;
    end

end
