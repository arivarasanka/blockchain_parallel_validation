%% heuristic_p1_simple.m
% Problem 1 simple transactions heuristic (height + volume + outdegree) on a ROW RANGE
clear; clc;

%% ---------- Config ----------
csvFile = 'tx_access_simple_21631019_21633879.csv';  % headered CSV
startId = 1;        % inclusive, 1-based row index
endId   = 20;     % inclusive, 1-based row index
p       = 2;        % capacity per round
% --------------------------------

%% ---------- Read needed columns & slice by range ----------
opts = detectImportOptions(csvFile, 'TextType','string');

needCols = ["tx_hash","access_read","access_write"];   
missing  = needCols(~ismember(needCols, string(opts.VariableNames)));
if ~isempty(missing)
    error('Missing required column(s): %s', strjoin(cellstr(missing), ', '));
end

% Force text columns to string (robust to empties/mixed)
for c = ["tx_hash","access_read","access_write"]
    idx = find(strcmp(opts.VariableNames, c), 1);
    if ~isempty(idx), opts.VariableTypes{idx} = 'string'; end
end
opts.SelectedVariableNames = intersect(opts.VariableNames, needCols);

Tfull = readtable(csvFile, opts);

% Clamp range to file bounds
Nfile  = height(Tfull);
startId = max(1, startId);
endId   = min(Nfile, endId);
if startId > endId
    error('Empty range: startId=%d > endId=%d (file has %d rows)', startId, endId, Nfile);
end

T = Tfull(startId:endId, :);
n = height(T);
fprintf('Row range [%d..%d] -> %d transactions selected.\n', startId, endId, n);
if n == 0, error('No transactions in the selected range.'); end

tx_hash      = string(T.tx_hash);          
access_read  = string(T.access_read);
access_write = string(T.access_write);

% Helper: normalize & tokenize access lists
parse_list = @(s) local_parse_list(s);

tStart=tic;

%% ---------- Build UNDIRECTED conflict graph from R/W sets ----------
% Conflict iff two txs share a key and at least one WRITES that key:
%   WW or WR or RW   (RR alone is NOT a conflict).
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

% writer-writer pairs + writer-reader pairs per key
wkeys = writers.keys;
for kk = 1:numel(wkeys)
    k = wkeys{kk};
    wlist = writers(k);

    % writer-writer (all unordered pairs)
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

    % writers vs readers
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
fprintf('Conflicts (undirected) on range: |E| = %d\n', size(E_und,1));

%% ---------- Convert conflicts to precedence arcs for P1 ----------
% Direct each {i,j} as i->j when i<j (row order) → acyclic precedence.
if isempty(E_und)
    E = zeros(0,2);
else
    E = E_und;   % already i<j by min/max above
end
m = size(E,1);
fprintf('Precedence arcs built: %d\n', m);

%% ---------- Build adjacency, indegree, outdegree ----------
succ   = cell(n,1);     % successors per node
pred   = cell(n,1);     % predecessors per node
indeg  = zeros(n,1);
outdeg = zeros(n,1);

for k = 1:m
    i = E(k,1); j = E(k,2);
    succ{i}(end+1) = j;        
    pred{j}(end+1) = i;        
    indeg(j)  = indeg(j)  + 1;
    outdeg(i) = outdeg(i) + 1;
end

%% ---------- Heights h(i) & Volumes v(i) ----------
% h(i) = 1 + max_{j in succ(i)} h(j); sinks have h=1
% v(i) = 1 + sum_{j in succ(i)} v(j); sinks have v=1  (double-counts across merges)
h = ones(n,1);
v = ones(n,1);
% Backward sweep okay because edges go forward in index (i<j).
for i = n:-1:1
    if ~isempty(succ{i})
        si = succ{i};
        h(i) = 1 + max(h(si));
        v(i) = 1 + sum(v(si));
    end
end

%% ---------- Heuristic scheduling ----------
% Key K(i) = (-h(i), -v(i), -outdeg(i), +id(i)); pick p smallest (lexicographic).
assignedRound = zeros(n,1);
remaining = n;

% Initialize ready set: indeg == 0
ready = find(indeg == 0).';

roundNum = 0;
while remaining > 0
    roundNum = roundNum + 1;

    % Order ready by key (ascending): [-h, -v, -outdeg, +id]
    if ~isempty(ready)
        K = [ -h(ready), -v(ready), -outdeg(ready), ready(:) ];
        [~,ord] = sortrows(K, [1 2 3 4]);   % lexicographic
        ready = ready(ord);
    end

    % Pick up to p jobs
    take = min(p, numel(ready));
    picked = ready(1:take);
    ready(1:take) = [];                 % remove picked from ready

    % Assign to this round
    for u = picked
        assignedRound(u) = roundNum;
    end
    remaining = remaining - numel(picked);

    % Update indegrees and ready set
    for u = picked
        Su = succ{u};
        for t = 1:numel(Su)
            vtx = Su(t);
            indeg(vtx) = indeg(vtx) - 1;
            if indeg(vtx) == 0
                ready(end+1) = vtx;     
            end
        end
    end

    % Safety (shouldn't trigger for a DAG)
    if numel(picked) == 0
        warning('No ready jobs found in round %d (graph may have a cycle).', roundNum);
        break;
    end
end

R_used = max(assignedRound);
fprintf('Heuristic finished: used rounds = %d\n', R_used);

%% ---------- Build xMat (n x R_used) & y ----------
xMat = zeros(n, R_used);
for i = 1:n
    r = assignedRound(i);
    if r >= 1
        xMat(i, r) = 1;
    end
end
y = any(xMat>0, 1).';
objval = sum(y);

%% ---------- Print schedule sample ----------
disp('--- Assignments (Transaction -> Round) ---');
Tshow = table( (1:n).', assignedRound, 'VariableNames', {'Transaction','Round'} );
disp( Tshow(1:min(n,40), :) );

loadPerRound = sum(xMat,1);
fprintf('Load per round (first min(20,R_used)): %s\n', mat2str(loadPerRound(1:min(20,R_used))));

totalSec = toc(tStart);
fprintf('Total time: %.3f s\n', totalSec);


%% ---------- Validation checks ----------
ok = true; tol = 1e-12;

% 1) exactly one per tx
rowSums = sum(xMat,2);
if any(abs(rowSums - 1) > tol)
    ok = false;
    bad = find(abs(rowSums - 1) > tol);
    fprintf('Exactly-one violated for %d tx (showing up to 10): %s\n', ...
        numel(bad), mat2str(bad(1:min(10,end))'));
end

% 2) capacity per round
over = find(loadPerRound > p + tol);
if ~isempty(over)
    ok = false;
    fprintf('Capacity overflow in rounds (showing up to 10): %s\n', ...
        mat2str(over(1:min(10,end))));
end

% 3) precedence: round(i) < round(j) for all (i->j)
viol_prec = [];
if ~isempty(E)
    ri = assignedRound(E(:,1));
    rj = assignedRound(E(:,2));
    bad = find(ri >= rj);
    if ~isempty(bad)
        ok = false;
        viol_prec = [E(bad,1), E(bad,2), ri(bad), rj(bad)];
        k = min(10, size(viol_prec,1));
        fprintf('Precedence violated (showing %d/%d): [i j round(i) round(j)]\n', k, size(viol_prec,1));
        disp( array2table(viol_prec(1:k,:), 'VariableNames', {'i','j','round_i','round_j'}) );
    end
end

% 4) objective consistency (no gaps ⇒ sum(y) == max round)
if objval ~= R_used
    ok = false;
    fprintf('Objective mismatch: sum(y)=%d vs used Rounds=%d\n', objval, R_used);
end

if ok
    fprintf('Verification: OK — one-per-tx, capacity OK, precedence OK, objective OK.\n');
end

%% ---------- Local parser ----------
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
