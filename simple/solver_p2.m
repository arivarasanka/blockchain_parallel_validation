%% solver_p2.m — P2 on a tx-id slice, conflicts from access_read/write
clear; clc;

%% -------- Config --------
csvFile = 'tx_access_simple_21631019_21635079.csv';
startId = 134001;          % 1-based row index (inclusive)
endId   = 140400;       % 1-based row index (inclusive)
p       = 8;           % capacity per round
R       = 100;         % number of rounds (e.g., floor(B/21000))
wCol    = 'weight_theoretical_eth';  % weight column to use
% --------------------------------

%% -------- Read needed columns robustly --------
opts = detectImportOptions(csvFile, 'TextType','string');
needCols = ["tx_hash","blockNumber","access_read","access_write",wCol];
have = ismember(needCols, string(opts.VariableNames));
if ~all(have)
    missing = needCols(~have);
    error('Missing required column(s) in CSV: %s', strjoin(missing, ', '));
end

% Force text columns to string
for c = ["tx_hash","access_read","access_write"]
    idx = find(strcmp(opts.VariableNames, c), 1);
    opts.VariableTypes{idx} = 'string';
end
% Select only the needed columns
opts.SelectedVariableNames = cellstr(needCols);

Tfull = readtable(csvFile, opts);
Nfull = height(Tfull);

% Bound and slice by row ids
startId = max(1, min(startId, Nfull));
endId   = max(1, min(endId,   Nfull));
if endId < startId
    error('endId (%d) < startId (%d).', endId, startId);
end
T = Tfull(startId:endId, :);
n = height(T);

fprintf('Slice [%d..%d] → %d transactions.\n', startId, endId, n);
if n == 0, error('Empty slice.'); end

tx_hash = string(T.tx_hash);  % may be empty strings; used only for display

%% -------- Weights (from wCol) --------
wraw = T.(wCol);
if isstring(wraw) || iscellstr(wraw)
    w = str2double(string(wraw));
else
    w = double(wraw);
end
w(~isfinite(w)) = 0;   % sanitize any NaNs/Infs
w = double(w(:));      % n-by-1

%% -------- Parse access lists & build conflict graph (undirected) --------
access_read  = string(T.access_read);
access_write = string(T.access_write);

parse_list = @(s) local_parse_list(s);

writers = containers.Map('KeyType','char','ValueType','any');  % key -> [local idx ...]
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

% Build undirected conflict edges: (WW) and (WR/RW). RR alone is allowed.
Ei = []; Ej = [];
wkeys = writers.keys;
for kk = 1:numel(wkeys)
    k = wkeys{kk};
    wlist = writers(k);                % local tx indices that WRITE key k

    % writer-writer pairs
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

    % writer-reader pairs
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

E = [Ei, Ej];
if ~isempty(E)
    E = unique(E, 'rows');  
end
fprintf('Conflicts built on slice: |E| = %d\n', size(E,1));

%% -------- Solve P2 --------
[xMat, objval, exitflag, output] = solver_min_rounds_p2(w, E, R, p);

usedRounds = nnz(sum(xMat,1));
fprintf('exitflag=%d | total weight=%.6g | used rounds=%d / %d | selected=%d\n', ...
        exitflag, objval, usedRounds, R, nnz(any(xMat,2)));

%% -------- Build in-memory schedule table (no file I/O) --------
rows = R * p;
R_col = zeros(rows,1);
P_col = zeros(rows,1);
IXcol = zeros(rows,1);      % local row index (1..n)
W_col = zeros(rows,1);
H_col = strings(rows,1);    % tx_hash (if available)

t = 0;
for r = 1:R
    idx = find(xMat(:,r) > 0.5);           % local indices in round r
    % Optional: sort by weight in-round
    % [~,ord] = sort(w(idx),'descend'); idx = idx(ord);
    for s = 1:p
        t = t + 1;
        R_col(t) = r;
        P_col(t) = s;
        if s <= numel(idx)
            i = idx(s);
            IXcol(t) = i;
            W_col(t) = w(i);
            H_col(t) = tx_hash(i);
        else
            IXcol(t) = 0;
            W_col(t)  = 0;
            H_col(t)  = "";
        end
    end
end

Sched = table(R_col, P_col, IXcol, H_col, W_col, ...
              'VariableNames', {'R','p','tx_local','tx_hash','W'});
disp(Sched(Sched.tx_local>0, :));

%% -------- Verification --------
tol = 1e-9;

% 1) At-most-one per tx
assignPerTx = sum(xMat, 2);
viol_tx = find(assignPerTx > 1 + tol);

% 2) Capacity per round
loadPerRound = sum(xMat, 1);
viol_cap = find(loadPerRound > p + tol);

% 3) No conflicts inside any round
badPairs = [];   % rows: [round, i, j]
if ~isempty(E)
    nE = size(E,1);
    for r = 1:R
        idx = find(xMat(:,r) > 0.5);
        if numel(idx) > 1
            sel = false(n,1); sel(idx) = true;
            hit = sel(E(:,1)) & sel(E(:,2));   % edges whose both ends are in round r
            if any(hit)
                eidx = find(hit);
                badPairs = [badPairs; [r*ones(numel(eidx),1), E(eidx,1), E(eidx,2)]];
            end
        end
    end
end

% 4) Objective consistency
selected = any(xMat > 0.5, 2);
obj_check = sum(w(selected));
obj_ok = abs(obj_check - objval) <= 1e-6 * max(1, abs(objval));

% 5) Report
if isempty(viol_tx) && isempty(viol_cap) && isempty(badPairs) && obj_ok
    fprintf('Verification: OK  | tx<=1 per job, capacity OK, no conflicts, objective matches.\n');
else
    fprintf('Verification: issues detected\n');
    if ~isempty(viol_tx)
        fprintf('  - At-most-one violated for %d tx (up to 10): %s\n', ...
            numel(viol_tx), mat2str(viol_tx(1:min(10,end))'));
    end
    if ~isempty(viol_cap)
        fprintf('  - Capacity violated in %d round(s) (up to 10): %s\n', ...
            numel(viol_cap), mat2str(viol_cap(1:min(10,end))));
    end
    if ~isempty(badPairs)
        k = min(10, size(badPairs,1));
        fprintf('  - Conflicting pairs inside rounds (showing %d/%d):\n', k, size(badPairs,1));
        disp(array2table(badPairs(1:k,:), 'VariableNames', {'Round','i','j'}));
    end
    if ~obj_ok
        fprintf('  - Objective mismatch: solver=%.6g vs recomputed=%.6g\n', objval, obj_check);
    end
end

%% -------- Local parser --------
function toks = local_parse_list(s)
    % Parse an access list cell into a vector of lowercase tokens.
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
