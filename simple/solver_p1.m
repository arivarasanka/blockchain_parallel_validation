%% Problem 1 (simple) — MILP Solver over a row-range [startId..endId]
clear; clc;

%% ---------- Config ----------
csvFile = 'tx_access_simple_21631019_21633879.csv';
startId = 68401;          % 1-based row index (inclusive)
endId   = 68800;       % 1-based row index (inclusive)
p       = 4;          % capacity per round
slackR  = 5;          % R = max(ceil(n/p), L2) + slack
% -----------------------------

%% ---------- Read needed columns ----------
opts = detectImportOptions(csvFile, 'TextType','string');
needCols = ["tx_hash","blockNumber","access_read","access_write"]; % blockNumber optional here
opts.SelectedVariableNames = intersect(opts.VariableNames, needCols);

% Force text cols to string (robust to empties)
for c = ["tx_hash","access_read","access_write"]
    if any(strcmp(opts.VariableNames, c))
        idx = find(strcmp(opts.VariableNames, c),1);
        opts.VariableTypes{idx} = 'string';
    end
end

Tfull = readtable(csvFile, opts);

% Clamp range to file bounds
Nfile = height(Tfull);
startId = max(1, startId);
endId   = min(Nfile, endId);
if startId > endId
    error('Empty range: startId=%d > endId=%d (file has %d rows).', startId, endId, Nfile);
end

% Slice window
T = Tfull(startId:endId, :);
n = height(T);
fprintf('Window [%d..%d] -> %d transactions\n', startId, endId, n);
if n == 0
    error('No transactions in the selected window.');
end

% Normalize to string columns for parsing
if ~ismember("access_read",  T.Properties.VariableNames),  error('Column "access_read" not found.');  end
if ~ismember("access_write", T.Properties.VariableNames),  error('Column "access_write" not found.'); end
access_read  = string(T.access_read);
access_write = string(T.access_write);

%% ---------- Helper: parse an access list string into tokens ----------
parse_list = @(s) local_parse_list(s);

%% ---------- Build conflict graph (undirected) ----------
% Conflict if share any key and at least one writes: W∩W, W∩R (R∩R allowed).
writers = containers.Map('KeyType','char','ValueType','any');
readers = containers.Map('KeyType','char','ValueType','any');

for i = 1:n
    Rset = parse_list(access_read(i));
    Wset = parse_list(access_write(i));

    % index writes
    for t = 1:numel(Wset)
        k = char(Wset(t));
        if isKey(writers,k), writers(k) = [writers(k), i];
        else,                writers(k) = i;
        end
    end
    % index reads
    for t = 1:numel(Rset)
        k = char(Rset(t));
        if isKey(readers,k), readers(k) = [readers(k), i];
        else,                readers(k) = i;
        end
    end
end

Ei = []; Ej = [];

% For each key: add writer-writer pairs and writer-reader pairs
wkeys = writers.keys;
for kk = 1:numel(wkeys)
    k = wkeys{kk};
    wlist = writers(k);             % tx indices that WRITE key k

    % writer-writer (all pairs)
    if numel(wlist) >= 2
        for a = 1:numel(wlist)-1
            ia = wlist(a);
            for b = a+1:numel(wlist)
                ib = wlist(b);
                i1 = min(ia,ib); i2 = max(ia,ib);
                Ei(end+1,1) = i1; Ej(end+1,1) = i2; 
            end
        end
    end

    % writers vs readers (cartesian)
    if isKey(readers, k)
        rlist = readers(k);
        for a = 1:numel(wlist)
            ia = wlist(a);
            for b = 1:numel(rlist)
                ib = rlist(b);
                if ia == ib, continue; end
                i1 = min(ia,ib); i2 = max(ia,ib);
                Ei(end+1,1) = i1; Ej(end+1,1) = i2; 
            end
        end
    end
end

E_und = [Ei, Ej];
if ~isempty(E_und), E_und = unique(E_und, 'rows'); end
fprintf('Conflicts (undirected) built: |E| = %d\n', size(E_und,1));

%% ---------- Convert conflicts to precedence arcs for P1 ----------
% P1 needs directed precedence consistent with the fixed order 1..n in THIS WINDOW.
% Direct each undirected conflict {i,j} as i->j when i<j (local indices).
if isempty(E_und)
    E = zeros(0,2);
else
    E = E_und;  % already i<j by min/max above
end

%% ---------- Compute R and solve P1 ----------
L1 = ceil(n / p);
L2 = longest_chain_length(n, E);
R  = max(L1, L2) + slackR;
fprintf('Lower bounds: ceil(n/p)=%d, longest_chain=%d → R=%d\n', L1, L2, R);

[xMat, y, objval, exitflag, output] = solve_min_rounds_p1(n, R, p, E); 

fprintf('Objective (sum used rounds) = %.0f\n', objval);
fprintf('exitflag = %d (1=optimal)\n', exitflag);

% Assigned round per tx (local id 1..n; original file row = startId + i - 1)
assignedRound = zeros(n,1);
for i = 1:n
    assignedRound(i) = find(xMat(i,:)>0.5, 1, 'first');
end

fprintf('--- First 40 assignments (LocalTx -> Round | FileRow) ---\n');
locIds  = (1:n).';
fileIds = startId - 1 + locIds;
show = table(locIds(1:min(40,n)), assignedRound(1:min(40,n)), fileIds(1:min(40,n)), ...
             'VariableNames', {'LocalTx','Round','FileRow'});
disp(show);

%% ---------- Verification ----------
tol = 1e-9;

% 1) exactly one per tx
rowSums = sum(xMat, 2);
viol_oneround = find(abs(rowSums - 1) > tol);

% 2) capacity & link to y
loadPerRound = sum(xMat, 1);
viol_capacity = find(loadPerRound > p + tol);
viol_link_on  = find((loadPerRound > tol) & (y.' < 1 - 1e-6));

% 3) no-gaps on y
if numel(y) >= 2
    viol_nogaps = find( y(1:end-1) - y(2:end) < -1e-9 );
else
    viol_nogaps = [];
end

% 4) precedence: i->j ⇒ round(i) < round(j)
viol_prec = [];
if ~isempty(E)
    lhs = assignedRound(E(:,1));
    rhs = assignedRound(E(:,2));
    bad = find(lhs >= rhs);
    if ~isempty(bad)
        viol_prec = [E(bad,1), E(bad,2), lhs(bad), rhs(bad)];
    end
end

% 5) objective consistency
lastUsedFromX = find(loadPerRound > tol, 1, 'last'); if isempty(lastUsedFromX), lastUsedFromX = 0; end
obj_from_y    = sum(y);
obj_ok = (abs(obj_from_y - objval) <= 1e-6*max(1,abs(objval))) && ...
         (obj_from_y == lastUsedFromX);

% Report
if isempty(viol_oneround) && isempty(viol_capacity) && isempty(viol_link_on) && ...
   isempty(viol_nogaps)  && isempty(viol_prec)      && obj_ok
    fprintf('Verification: OK  | one-per-tx, capacity OK, y-link OK, no-gaps OK, precedence OK, objective matches.\n');
else
    fprintf('Verification: issues detected\n');
    if ~isempty(viol_oneround)
        fprintf('  - Exactly-one violated (up to 10): %s\n', mat2str(viol_oneround(1:min(10,end))'));
    end
    if ~isempty(viol_capacity)
        fprintf('  - Capacity overflow rounds (up to 10): %s\n', mat2str(viol_capacity(1:min(10,end))));
    end
    if ~isempty(viol_link_on)
        fprintf('  - x>0 but y=0 in rounds (up to 10): %s\n', mat2str(viol_link_on(1:min(10,end))));
    end
    if ~isempty(viol_nogaps)
        fprintf('  - No-gaps violated at r indices (up to 10): %s\n', mat2str(viol_nogaps(1:min(10,end))));
    end
    if ~isempty(viol_prec)
        k = min(10, size(viol_prec,1));
        fprintf('  - Precedence violated (showing %d/%d): [i j round(i) round(j)]\n', k, size(viol_prec,1));
        disp( array2table(viol_prec(1:k,:), 'VariableNames', {'i','j','round_i','round_j'}) );
    end
    if ~obj_ok
        fprintf('  - Objective mismatch: sum(y)=%d, lastUsedFromX=%d, solver obj=%.0f\n', ...
            obj_from_y, lastUsedFromX, objval);
    end
end

%% ---------- Local function: parse list ----------
function toks = local_parse_list(s)
    % Robustly parse access list cell into a vector of lowercased tokens.
    if ismissing(s) || strlength(s) == 0
        toks = strings(0,1); return;
    end
    t = lower(strtrim(s));
    % strip common wrappers: brackets/braces/quotes
    t = replace(t, ["[","]","{","}","'",char(34)], "");
    % normalize separators to single spaces
    t = regexprep(t, '[;,\|\s]+', ' ');
    if strlength(t) == 0
        toks = strings(0,1); return;
    end
    toks = string(strsplit(t, ' '));
    toks = toks(toks ~= "");
end
