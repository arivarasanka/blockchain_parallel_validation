function [xMat, y, objval, exitflag, output] = solve_min_rounds_p1(n, R, p, E)
% Inputs: n (#transactions), R (max rounds), p (capacity), E (m-by-2 precedence arcs [i j])

% Decision vector z stacks all variables: z = [x(:); y]
numX = n*R;          % number of x(i,r) variables
numY = R;            % number of y(r) variables
Nvar = numX + numY;  % total variables.

% Helpers to map (i,r) and (r) into positions of z
idxX = @(i,r) (r-1)*n + i;   % x(i,r) -> 1..numX
idxY = @(r)   numX + r;      % y(r)   -> numX+1 .. numX+R

% Objective coefficients: minimize sum_r y_r
f = zeros(Nvar,1);
for r = 1:R
    f(idxY(r)) = 1;
end

% Allocate empty constraint matrices 
A = []; b = [];
Aeq = []; beq = [];

% Variable bounds (0 <= x,y <= 1)
lb = zeros(Nvar,1);
ub = ones(Nvar,1);

% STEP 3: Each transaction is assigned to exactly one round: sum_r x(i,r) = 1
Aeq_oneround = zeros(n, Nvar);   % n equalities, one per transaction i
beq_oneround = ones(n,1);        % right-hand side = 1 for each i

for i = 1:n
    for r = 1:R
        Aeq_oneround(i, idxX(i,r)) = 1;  % coefficient of x(i,r) in row i
    end
end

% Append to the global equality system
Aeq = [Aeq; Aeq_oneround];
beq = [beq; beq_oneround];

% STEP 4: Capacity-usage link per round
ii = zeros(n*R + R,1); jj = zeros(n*R + R,1); vv = zeros(n*R + R,1);
k = 0;
for r = 1:R
    for i = 1:n
        k=k+1; ii(k)=r; jj(k)=idxX(i,r); vv(k)=1;
    end
    k=k+1; ii(k)=r; jj(k)=idxY(r); vv(k)=-p;
end
A_caplink = sparse(ii(1:k), jj(1:k), vv(1:k), R, Nvar);
b_caplink = zeros(R,1);
A = [A; A_caplink]; b = [b; b_caplink];

% STEP 5: Precedence: x(i,r) + x(j,rp) <= 1  for all (i,j) in E, and all rp <= r
% For each (i,j) in E and r=1..R:
%   sum_{t=1}^r x(j,t) - sum_{t=1}^{r-1} x(i,t) <= 0
m = size(E,1);
ii = zeros(m*R*R,1); jj = zeros(m*R*R,1); vv = zeros(m*R*R,1);
row = 0; k = 0;
for e = 1:m
    i = E(e,1); j = E(e,2);
    for r_ = 1:R
        row = row + 1;
        for t = 1:r_       % + x(j,t)
            k=k+1; ii(k)=row; jj(k)=idxX(j,t); vv(k)=vv(k)+1;
        end
        for t = 1:(r_-1)   % - x(i,t)
            k=k+1; ii(k)=row; jj(k)=idxX(i,t); vv(k)=vv(k)-1;
        end
    end
end
A_prec = sparse(ii(1:k), jj(1:k), vv(1:k), m*R, Nvar);
b_prec = zeros(m*R,1);
A = [A; A_prec]; b = [b; b_prec];

% No-gaps constraint
A_nogaps = zeros(R-1, Nvar); b_nogaps = zeros(R-1,1);
for r = 1:R-1
    A_nogaps(r, idxY(r+1)) =  1;
    A_nogaps(r, idxY(r))   = -1;
end
A = [A; A_nogaps]; b = [b; b_nogaps];

% STEP 7: Declare integrality
intcon = 1:(numX + numY);   % x and y are integer (binary)


opts = optimoptions('intlinprog', ...
    'Display','iter', ...
    'Heuristics','advanced', ...
    'CutGeneration','intermediate');
tSolve = tic;

[z, objval, exitflag, output] = intlinprog( ...
    f, intcon, A, b, Aeq, beq, lb, ub, opts);

output.time_solver = toc(tSolve);   % << store solver wall time (seconds)


% Unpack solution into outputs
x = z(1:numX);
y = z(numX+1:end);

xMat = reshape(x, [n, R]);  % n-by-R assignment matrix

% STEP 8: Post-process (optional helpers)

% (a) make x explicitly 0/1 for neatness (handles tiny numeric tolerances)
xMat = round(xMat);

% (b) which rounds are used?
roundsUsed = find(y > 1e-6);

% (c) assigned round for each transaction i (since sum_r x(i,r) = 1)
assignedRound = zeros(n,1);
for i = 1:n
    rr = find(xMat(i,:) > 0.5, 1, 'first');
    if isempty(rr), rr = NaN; end
    assignedRound(i) = rr;
end
% assignedRound is not returned, but handy to inspect during debugging
% T = table((1:n).', assignedRound, 'VariableNames', {'Transaction','Round'});
% disp(T);

end
