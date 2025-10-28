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

function [xMat, s, e, M, objval, exitflag, output] = solve_min_rounds_p1(n, p, t, E, varargin)
% SOLVE_MIN_ROUNDS_P1  Heterogeneous, ordered-block scheduling MILP
% Minimize makespan M on p identical cores with:
% - Block-order precedence on conflicting txs:   s_j >= e_i for (i,j) in E
% - Assignment to exactly one core:              sum_c x(i,c) = 1
% - NEW: No overlap ON THE SAME CORE via z_ij gating:
%        w_ijc links to x, z_ij = sum_c w_ijc, then
%        e_i <= s_j + B[(1-y_ij) + (1-z_ij)]
%        e_j <= s_i + B[   y_ij  + (1-z_ij)]

% INPUTS
%   n         : #transactions (jobs) 1..n
%   p         : #cores
%   t         : n-by-1 vector of processing times t_i
%   E         : m-by-2 directed edges (i,j): conflicts with i<j in block
%
% OUTPUTS
%   xMat     : n-by-p binary assignment matrix (row i → chosen core)
%   s, e     : n-by-1 start/end times
%   M        : scalar makespan
%   objval   : optimal objective value (=M)
%   exitflag : intlinprog exit flag
%   output   : solver output struct

% ------------ Parse options ------------
ip = inputParser;
ip.addParameter('BigM',[]);
ip.addParameter('Display','iter');
ip.parse(varargin{:});
B = ip.Results.BigM;
dispLevel = ip.Results.Display;

t = t(:);
if isempty(B), B = sum(t); end

% ------------ Index maps for variables ------------
% Decision vector z = [ x(:) ; s ; e ; M ; y(:) ; w(:) ; zz(:) ]
numX = n*p;            % x(i,c)
numS = n;              % s(i)
numE = n;              % e(i)
numM = 1;              % M

idxX = @(i,c) (c-1)*n + i;             % x(i,c) in 1..numX
idxS = @(i)  numX + i;                 % s(i)
idxE = @(i)  numX + numS + i;          % e(i)
idxM =        numX + numS + numE + 1;  % M

% ------------ Build unordered pair list for non-conflicts ℘_non ------------
% adj(i,j)=true iff (i,j)∈E (directed). We want pairs with NO edge either way.
adj = false(n,n);
if ~isempty(E)
    adj(sub2ind([n,n], E(:,1), E(:,2))) = true;
end

pairsI = []; pairsJ = [];
for i = 1:n-1
    js = (i+1):n;
    % keep only non-edges in BOTH directions
    mask = ~(adj(i,js) | adj(js,i).');   % 1×|js| logical row
    js = js(mask);
    if ~isempty(js)
        pairsI = [pairsI; i*ones(numel(js),1)];
        pairsJ = [pairsJ; js(:)];
    end
end
mPairs = numel(pairsI);

% (Optional) sanity check
% fprintf('Using %d non-edge pairs for z-gated constraints.\n', mPairs);


% y_ij for each pair in ℘_non
idxY = @(k) idxM + k;                 % k=1..mPairs
lastY = idxY(mPairs);

% w_ijc for each pair and core
% Linearize index: W(k,c) laid out with c fastest
idxW = @(k,c) lastY + (k-1)*p + c;    % k=1..mPairs, c=1..p
lastW = idxW(mPairs, p);

% z_ij for each pair
idxZ = @(k) lastW + k;                % k=1..mPairs
Nvar  = idxZ(mPairs);                 % total vars

% ------------ Objective: minimize M ------------
f = zeros(Nvar,1);
f(idxM) = 1;

% ------------ Bounds & integrality ------------
lb = zeros(Nvar,1);  % s,e,M ≥ 0; binaries ≥ 0
ub = ones(Nvar,1);   % default upper bound 1
% s,e,M continuous and unbounded above:
ub(idxS(1):idxS(n)) = inf;
ub(idxE(1):idxE(n)) = inf;
ub(idxM) = inf;

% Integer vars: x, y, w, z
intcon = [ (1:numX), (idxY(1):idxY(mPairs)), (lastY+1:lastW), (idxZ(1):idxZ(mPairs)) ];

% ------------ Equality constraints ------------
% (1) sum_c x(i,c) = 1
Aeq1 = sparse(n, Nvar);
beq1 = ones(n,1);
for i = 1:n
    for c = 1:p
        Aeq1(i, idxX(i,c)) = 1;
    end
end

% (2) e_i - s_i = t_i
Aeq2 = sparse(n, Nvar);
beq2 = t(:);
for i = 1:n
    Aeq2(i, idxE(i)) =  1;
    Aeq2(i, idxS(i)) = -1;
end

% (3) z_ij = sum_c w_ijc  (for all pairs)
Aeq3 = sparse(mPairs, Nvar);
beq3 = zeros(mPairs,1);
for k = 1:mPairs
    Aeq3(k, idxZ(k)) = 1;
    for c = 1:p
        Aeq3(k, idxW(k,c)) = -1;
    end
end

Aeq = [Aeq1; Aeq2; Aeq3];
beq = [beq1;  beq2;  beq3];

% ------------ Inequality constraints ------------
rowsI = {}; rowsJ = {}; rowsV = {}; rhs = [];
rowCounter = 0;

% helper: append one row sum a_j z_j <= b
    function appendRow(cols, vals, bval)
        rowCounter = rowCounter + 1;
        rowsI{rowCounter} = 1:numel(cols);
        rowsJ{rowCounter} = cols(:).';
        rowsV{rowCounter} = vals(:).';
        rhs(rowCounter,1) = bval;
    end

% (A) Precedence for directed edges: e_i - s_j <= 0
if ~isempty(E)
    for r = 1:size(E,1)
        i = E(r,1); j = E(r,2);
        appendRow([idxE(i), idxS(j)], [1, -1], 0);
    end
end

% (B) Makespan: e_i - M <= 0
for i = 1:n
    appendRow([idxE(i), idxM], [1, -1], 0);
end

% (C) Link w_ijc to x (for all non-conflict pairs and cores)
%     w_ijc <= x_ic,  w_ijc <= x_jc,  w_ijc >= x_ic + x_jc - 1  (⇒ -x_ic - x_jc + w_ijc <= -1)
% (C) Link w_ijc to x (for all non-conflict pairs and cores)
for k = 1:mPairs
    i = pairsI(k); j = pairsJ(k);
    for c = 1:p
        wkc = idxW(k,c); xic = idxX(i,c); xjc = idxX(j,c);

        % w_ijc <= x_ic
        appendRow([wkc, xic], [1, -1], 0);

        % w_ijc <= x_jc
        appendRow([wkc, xjc], [1, -1], 0);

        % w_ijc >= x_ic + x_jc - 1  ⇔  x_ic + x_jc - w_ijc <= 1   (CORRECTED)
        appendRow([xic, xjc, wkc], [1, 1, -1], 1);
    end
end


% (D) Gated same-core disjunction (NEW):
%     e_i <= s_j + B[(1-y_ij) + (1-z_ij)]  ⇒  e_i - s_j + B*y_ij + B*z_ij <= 2B
%     e_j <= s_i + B[   y_ij  + (1-z_ij)]  ⇒  e_j - s_i - B*y_ij + B*z_ij <=  B
for k = 1:mPairs
    i = pairsI(k); j = pairsJ(k);
    yk = idxY(k); zk = idxZ(k);
    appendRow([idxE(i), idxS(j), yk, zk], [1, -1,  B,  B], 2*B);
    appendRow([idxE(j), idxS(i), yk, zk], [1, -1, -B,  B],   B);
end

% ---- Assemble sparse A, b
if isempty(rowsI)
    A = sparse(0, Nvar); b = zeros(0,1);
else
    totalRows = numel(rhs);
    nnzTot = 0;
    for r = 1:totalRows, nnzTot = nnzTot + numel(rowsI{r}); end
    II = zeros(nnzTot,1); JJ = zeros(nnzTot,1); VV = zeros(nnzTot,1);
    pos = 0;
    for r = 1:totalRows
        m = numel(rowsI{r});
        II(pos+1:pos+m) = r;
        JJ(pos+1:pos+m) = rowsJ{r};
        VV(pos+1:pos+m) = rowsV{r};
        pos = pos + m;
    end
    A = sparse(II, JJ, VV, totalRows, Nvar);
    b = rhs;
end

% ------------ Solve ------------
opts = optimoptions('intlinprog','MaxTime',6000, 'Display', dispLevel);
[z, objval, exitflag, output] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts);

% ------------ Unpack solution ------------
x = z(1:numX);
s = z(idxS(1):idxS(n));
e = z(idxE(1):idxE(n));
M = z(idxM);

xMat = reshape(x, [n, p]);
xMat = round(xMat);  % cleanup tiny tolerances
end
