function L = longest_chain_length(n,E)
% returns length of the longest path in DAG (nodes 1..n)
if isempty(E), L = 1; return; end
G = cell(n,1); indeg = zeros(n,1);
for k=1:size(E,1)
    i=E(k,1); j=E(k,2);
    G{i}(end+1)=j; indeg(j)=indeg(j)+1;
end
q = find(indeg==0).';  % topo queue
dist = ones(n,1);      % path length ending at node (at least 1)
h = 1;
while h <= numel(q)
    u = q(h); h = h+1;
    for v = G{u}
        dist(v) = max(dist(v), dist(u)+1);
        indeg(v) = indeg(v)-1;
        if indeg(v)==0, q(end+1)=v; end
    end
end
L = max(dist);
end
