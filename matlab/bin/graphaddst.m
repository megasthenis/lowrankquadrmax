function graph = graphaddst(G, fromS, toT)
% G : the adjacency matrix of a graph on n vertices.


[n, j] = size(G);

if n ~= j 
   error('G not square'); 
end

S = n+1;
T = n+2;

if nargin < 2
    fromS = 1:n;
    toT = 1:n;
end

assert(~any(fromS > n) && ~any(fromS < 1))
assert(~any(toT > n) && ~any(toT < 1))

graph.adj = sparse(n+2, n+2);
graph.adj(1:n, 1:n) = G;
graph.adj(S, fromS) = 1;
graph.adj(toT, T) = 1;

graph.n = n + 2;            % Number of vertices (eacy access)
graph.m = nnz(graph.adj);   % Number of arcs (eacy access)
graph.s = S;                % i.d. of node s
graph.t = T;                % i.d. of node t
