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

function [xMat, objval, exitflag, output] = solver_min_rounds_p2(w, E, R, p)
% SOLVER_MIN_ROUNDS_P2  Maximum-weight selection with conflicts, R rounds, cap p per round
% Inputs:
%   w : n-by-1 weights (w_i >= 0)
%   E : m-by-2 int matrix of conflicts (1..n), undirected pairs
%   R : number of rounds (colors)
%   p : capacity per round (max #tx per round)
% Outputs:
%   xMat     : n-by-R binary assignment (x(i,r)=1 if tx i in round r)
%   objval   : optimal total weight
%   exitflag : intlinprog status (1 = optimal)
%   output   : solver diagnostics

w = double(w(:));
n = numel(w);
if isempty(E), E = zeros(0,2); end
E = double(E);

% ------------- Variables and indexing -------------
numX = n*R;     % only x(i,r) binaries
Nvar = numX;
idxX = @(i,r) (r-1)*n + i;

% ------------- Objective: maximize sum_i,r w_i x(i,r) -------------
% intlinprog minimizes, so use -w block-repeated across rounds
f = sparse(Nvar,1);
for r = 1:R
    f((r-1)*n + (1:n)) = -w;
end

% ------------- Bounds & integrality -------------
lb = zeros(Nvar,1);
ub = ones(Nvar,1);
intcon = 1:Nvar;

% ------------- Inequalities A*z <= b -------------
% (1) At-most-one per transaction: sum_r x(i,r) <= 1   [n rows]
ii = zeros(n*R,1); jj = zeros(n*R,1); vv = ones(n*R,1);
k = 0;
for r = 1:R
    base = (r-1)*n;
    idx = k+1 : k+n;
    ii(idx) = (1:n).';
    jj(idx) = base + (1:n).';
    k = k + n;
end
A_once = sparse(ii, jj, vv, n, Nvar);
b_once = ones(n,1);

% (2) Capacity per round: sum_i x(i,r) <= p              [R rows]
ii = zeros(n*R,1); jj = zeros(n*R,1); vv = ones(n*R,1);
k = 0;
for r = 1:R
    base = (r-1)*n;
    idx = k+1 : k+n;
    ii(idx) = r;
    jj(idx) = base + (1:n).';
    k = k + n;
end
A_cap = sparse(ii, jj, vv, R, Nvar);
b_cap = p*ones(R,1);

% (3) Same-round conflicts: x(i,r) + x(j,r) <= 1         [|E|*R rows]
if isempty(E)
    A_conf = sparse(0, Nvar); b_conf = zeros(0,1);
else
    m = size(E,1);
    ii = zeros(2*m*R,1); jj = zeros(2*m*R,1); vv = ones(2*m*R,1);
    row = 0; k = 0;
    for e = 1:m
        i = E(e,1); j = E(e,2);
        for r = 1:R
            row = row + 1;
            k = k + 1; ii(k) = row; jj(k) = idxX(i,r);
            k = k + 1; ii(k) = row; jj(k) = idxX(j,r);
        end
    end
    A_conf = sparse(ii(1:k), jj(1:k), vv(1:k), row, Nvar);
    b_conf = ones(row,1);
end

% Assemble
A = [A_once; A_cap; A_conf];
b = [b_once; b_cap; b_conf];

% ------------- Equalities (none) -------------
Aeq = sparse(0, Nvar);
beq = zeros(0,1);

% ------------- Solve -------------
opts = optimoptions('intlinprog', ...
    'Display','iter', ...
    'MaxTime', inf, ...             
    'RelativeGapTolerance', 1e-3);    

[z, fval, exitflag, output] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, opts);

if isempty(z)
    error('No solution returned: %s (exitflag=%d)', output.message, exitflag);
end

% ------------- Unpack solution -------------
x = z(1:numX);
xMat = reshape(x, [n, R]);    % i by r
xMat = round(xMat);           % clean tolerances

objval = -fval;               % back to max-weight
end
