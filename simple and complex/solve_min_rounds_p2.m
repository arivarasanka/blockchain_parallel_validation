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

function [xMat, s, e, objval, exitflag, output] = solve_min_rounds_p2(n, p, t, E, w, varargin)

ip = inputParser;
ip.addParameter('R', []);
ip.addParameter('Display', 'iter');
ip.parse(varargin{:});
B = ip.Results.R;
dispLevel = ip.Results.Display;

t = t(:); w = w(:);
if isempty(B), B = sum(t); end
assert(numel(t)==n && numel(w)==n, 't and w must be n-vectors');

% Upper-tri pairs (i<j)
pairsI = []; pairsJ = [];
for i = 1:n-1
    js = (i+1):n;
    pairsI = [pairsI; i*ones(numel(js),1)];
    pairsJ = [pairsJ; js(:)];
end
mPairs = numel(pairsI);

% Conflict mask (upper triangle)
isConf = false(n,n);
if ~isempty(E)
    for r = 1:size(E,1)
        i = E(r,1); j = E(r,2);
        if i==j, continue; end
        if i<j, isConf(i,j) = true; else, isConf(j,i) = true; end
    end
end

% Variable blocks
numX   = n*p;          % x(i,c)
numS   = n;            % s(i)
numE   = n;            % e(i)
numV   = n;            % v(i)
numY   = mPairs;       % y(i,j)
numZpr = mPairs;       % z_ij
numW   = mPairs*p;     % w_ijc

idxX  = @(i,c)                (c-1)*n + i;
offS  = numX;
idxS  = @(i)                  offS + i;
offE  = offS + numS;
idxE  = @(i)                  offE + i;
offV  = offE + numE;
idxV  = @(i)                  offV + i;
offY  = offV + numV;
idxY  = @(k)                  offY + k;
offZp = offY + numY;
idxZp = @(k)                  offZp + k;
offW  = offZp + numZpr;
idxW  = @(k,c)                offW + (k-1)*p + c;
Nvar  = offW + numW;

% Objective: maximize sum_i w_i v_i  (minimize -w)
f = zeros(Nvar,1);
for i = 1:n, f(idxV(i)) = -w(i); end

% Bounds
lb = zeros(Nvar,1);
ub = ones(Nvar,1);
ub(idxS(1):idxS(n)) = inf;   % s continuous
ub(idxE(1):idxE(n)) = inf;   % e continuous

% Integer variables (FIX HERE: contiguous W block)
intcon = [ 1:numX, ...                             % x
           idxV(1):idxV(n), ...                   % v
           idxY(1):idxY(mPairs), ...              % y
           idxZp(1):idxZp(mPairs), ...            % z_ij
           (offW+1):Nvar ];                       % all w_ijc

% Equalities
% (A) sum_c x_ic = v_i
AeqA = spalloc(n, Nvar, n*p + n);
beqA = zeros(n,1);
for i = 1:n
    for c = 1:p, AeqA(i, idxX(i,c)) = 1; end
    AeqA(i, idxV(i)) = -1;
end

% (B) e_i - s_i = t_i v_i
AeqB = spalloc(n, Nvar, 3*n);
beqB = zeros(n,1);
for i = 1:n
    AeqB(i, idxE(i)) =  1;
    AeqB(i, idxS(i)) = -1;
    AeqB(i, idxV(i)) = -t(i);
end

% (C) z_ij = sum_c w_ijc  (only for non-conflicts)
numEqC = sum(~isConf(sub2ind([n,n], pairsI, pairsJ)));
AeqC = spalloc(numEqC, Nvar, numEqC*(p+1));
beqC = zeros(numEqC,1);
rowC = 0;
for k = 1:mPairs
    i = pairsI(k); j = pairsJ(k);
    if isConf(i,j), continue; end
    rowC = rowC + 1;
    AeqC(rowC, idxZp(k)) = 1;
    for c = 1:p
        AeqC(rowC, idxW(k,c)) = -1;
    end
end

Aeq = [AeqA; AeqB; AeqC];
beq = [beqA;  beqB;  beqC];

% Inequalities builder
rowsII = {}; rowsJJ = {}; rowsVV = {}; rhs = [];
rowCounter = 0;
    function addRow(cols, vals, bval)
        rowCounter = rowCounter + 1;
        rowsII{rowCounter} = ones(1,numel(cols))*rowCounter;
        rowsJJ{rowCounter} = cols(:).';
        rowsVV{rowCounter} = vals(:).';
        rhs(rowCounter,1) = bval;
    end

% Activation upper bounds: s_i <= B v_i ; e_i <= B v_i
for i = 1:n
    addRow([idxS(i), idxV(i)], [1, -B], 0);
    addRow([idxE(i), idxV(i)], [1, -B], 0);
end

% Conflict constraints (i,j) in E, i<j
for k = 1:mPairs
    i = pairsI(k); j = pairsJ(k);
    if ~isConf(i,j), continue; end
    yk = idxY(k);
    % e_i - s_j <= B(1 - y_ij)  -> e_i - s_j + B*y_ij <= B
    addRow([idxE(i), idxS(j), yk], [1, -1,  B],  B);
    % e_j - s_i <= B*y_ij       -> e_j - s_i - B*y_ij <= 0
    addRow([idxE(j), idxS(i), yk], [1, -1, -B],  0);
end

% Non-conflict gating: w_ijc <= x_ic; w_ijc <= x_jc; -w_ijc + x_ic + x_jc <= 1
for k = 1:mPairs
    i = pairsI(k); j = pairsJ(k);
    if isConf(i,j), continue; end
    for c = 1:p
        addRow([idxW(k,c), idxX(i,c)], [1, -1], 0);
        addRow([idxW(k,c), idxX(j,c)], [1, -1], 0);
        addRow([idxW(k,c), idxX(i,c), idxX(j,c)], [-1, 1, 1], 1);
    end
end

% Non-conflict ordering:
% e_i - s_j + B*y_ij + B*z_ij <= 2B
% e_j - s_i - B*y_ij + B*z_ij <= B
for k = 1:mPairs
    i = pairsI(k); j = pairsJ(k);
    if isConf(i,j), continue; end
    yk = idxY(k); zk = idxZp(k);
    addRow([idxE(i), idxS(j), yk, zk], [1, -1,  B,  B], 2*B);
    addRow([idxE(j), idxS(i), yk, zk], [1, -1, -B,  B],   B);
end

% Assemble A,b
if rowCounter == 0
    A = sparse(0, Nvar); b = zeros(0,1);
else
    nnzTot = 0;
    for r = 1:rowCounter, nnzTot = nnzTot + numel(rowsJJ{r}); end
    II = zeros(nnzTot,1); JJ = zeros(nnzTot,1); VV = zeros(nnzTot,1);
    pos = 0;
    for r = 1:rowCounter
        m = numel(rowsJJ{r});
        II(pos+1:pos+m) = rowsII{r};
        JJ(pos+1:pos+m) = rowsJJ{r};
        VV(pos+1:pos+m) = rowsVV{r};
        pos = pos + m;
    end
    A = sparse(II, JJ, VV, rowCounter, Nvar);
    b = rhs;
end

% Solve
opts = optimoptions('intlinprog', 'Display', dispLevel, 'MaxTime', 6000);
[z, objval, exitflag, output] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts);

% Unpack
objval = -objval;
x = z(1:numX);
s = z(idxS(1):idxS(n));
e = z(idxE(1):idxE(n));
xMat = reshape(round(x), [n, p]);
end
