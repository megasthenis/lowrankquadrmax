function [x] = rank1solver_stpath(w, options)
% Structered PCA - st path on a graph: Solve the rank-1 problem. 
% Returns unit norm vector x of dimension equal to w.
% The support of x must coincide with vertices lying accross an st path
% of a graph G specified in the options.

%global TOL;

n = numel(w);
G = options.graph;

% in-degree of each vertex in G
indegree = sum( (G.adj ~= 0), 1);

% arc (u, v) [u,v not in s,t] gets the weight of its destination vertex v.
% arc (s, u) gets weight 0 for all u.
% arc (u, t) gets weight 0 for all u.
weights = zeros(G.m, 1);    % arc weights
ctr = 0;
for v = 1:G.n    
    if v == G.s || v == G.t
        val = 0;
    else
        val = w(v)^2;
    end
    weights(ctr+1:ctr+indegree(v)) = val;
    ctr = ctr + indegree(v);
end

%assert(numel(edge_weights) == nnz(G.adj))
[~, path, ~] = graphshortestpath(G.adj, G.s, G.t,...
                                'Directed', true, ...
                                'Method', 'Acyclic', ...
                                'Weights', -weights );
assert(path(1)   == G.s)
assert(path(end) == G.t)

support = path(2:end-1);    % remove s and t from optimal path

x = sparse(support, 1, ...  % support (subvector indices)
           w(support), ...  % nonzero values of x
           n, 1 ...         % size of x
    );

x = x / norm(x);

end % end of function
