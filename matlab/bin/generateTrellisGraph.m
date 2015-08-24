function [G, layerSizes] = generateTrellisGraph(parameters)
% Generate a directed acyclic graph on n vertices with the following 
% structure: 
% The n vertices are partitioned into L layers (labeled 1,...L),
% Each vertex of layer i has outgoing arcs towards all vertices of layer 


% FUNCTION CONSTANTS
FULL    = 0;
RANDOM  = 1;

%% Parse Input Parameters
if ~isfield(parameters, 'n')
    error('n: number of vertices not specified');
end

n = parameters.n;       % number of vertices

if n < 0 || floor(n) ~= n
    error('n: number of vertices must be positive integer');
end

if ~isfield(parameters, 'L')
    error('L: number of layers not specified');
end

if numel(parameters.L) == 1
    % user determined number of layers
    numLayers  = parameters.L;
    layerSizes = [];
elseif numel(parameters.L) > 1
    % user determined number of nodes per layer
    % Total number of vertices in layers must much total number of vertices
    assert(sum(parameters.L) == n)
    numLayers  = numel(parameters.L);
    layerSizes = parameters.L;
end


% if user does not specify out_degree, assume full trellis
if ~isfield(parameters, 'out_degree')
    mode = FULL;
else
    degree = parameters.out_degree;
    if degree < 0 || floor(degree) ~= degree
        warning('out_degree must be positive integer');
    end
    mode = RANDOM;
end


%% Layers

% Determine the size of each of the L layers
if isempty(layerSizes)
    k = ceil(n / numLayers);
    layerSizes = zeros(numLayers,1);
    remaining_vertices = n;
    for i = 1:numLayers
        layerSizes(i) = min([k, remaining_vertices]);
        remaining_vertices = remaining_vertices - layerSizes(i);
        if remaining_vertices == 0
           break; 
        end
    end
    if i < numLayers
       warning('Created graph with fewer layers.');
       numLayers = i;
    end
    layerSizes = layerSizes(1:numLayers);
end

%% Edges

% Determine the maximum number of edges
switch(mode)
case FULL
    m = sum(layerSizes(1:end-1) .* layerSizes(2:end));
case RANDOM
    m = sum(degree * layerSizes(1:end-1));
end

edges = zeros(m, 2); % from-to pairs

switch(mode)
case FULL
    
    % From each vertex to all vertices in next layer
    i = 1; offset = 0;
    for l = 1:numLayers-1
        currlayer = offset+1 : offset+layerSizes(l);
        nextlayer = offset+layerSizes(l)+1 : offset+layerSizes(l)+layerSizes(l+1);
        deg       = numel(nextlayer);
        for from = currlayer
           edges(i: i+deg-1, 1) = from;
           edges(i: i+deg-1, 2) = nextlayer;
           i = i + deg;
        end
        offset = offset + layerSizes(l);
    end
    
case RANDOM
    
    % From each vertex to "degree" randomly chosen vertices in next layer
    % (if degree nodes exist)
    i = 1; offset = 0;
    for l = 1:numLayers-1
        currlayer = offset+1 : offset+layerSizes(l);
        nextlayer = offset+layerSizes(l)+1 : offset+layerSizes(l)+layerSizes(l+1);
        cand = numel(nextlayer);
        deg = min([degree, cand]);
        
        for from = currlayer
           ro = randperm(cand);
           edges(i: i+deg-1, 1) = from;
           edges(i: i+deg-1, 2) = nextlayer(ro(1:deg));
           i = i + deg;
        end
        offset = offset + layerSizes(l);
    end
    
end



G = sparse(edges(1:i-1, 1), edges(1:i-1, 2),... % edge positions
            ones(i-1,1), ...                      % edge value
            n, n);                              % adjacency size



end % end of function


