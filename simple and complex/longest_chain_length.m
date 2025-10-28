% -------------------------------------------------------------------------
% Author: [Arivarasan Karmegam]
% Supervisor Name: [Antonio Fernandez Anta]
% Copyright (C) [2025] [IMDEA Networks Institute]. All rights reserved.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% -------------------------------------------------------------------------

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
