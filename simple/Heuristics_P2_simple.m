%% heuristic_p2_simple.m — Problem 2 heuristic (score = w/(deg+1), greedy packing)
clear; clc;

%% -------- Inputs --------
csvFile = 'tx_access_simple_21631019_21635079.csv';
startId = 90801;          % 1-based row index (inclusive)
endId   = 92400;       % 1-based row index (inclusive)
p       = 2;         % capacity per round
R       = 100;         % number of rounds (e.g., floor(B/21000))
wCol    = 'weight_theoretical_eth';  % weight column to use
% -------------------------

%% -------- Read needed columns robustly --------
needCols = ["tx_hash","access_read","access_write", wCol];
opts = detectImportOptions(csvFile, 'TextType','string');
% Force present text cols to string for safety
for c = ["tx_hash","access_read","access_write"]
    if any(strcmp(opts.VariableNames, c))
        idx = find(strcmp(opts.VariableNames, c),1);
        opts.VariableTypes{idx} = 'string';
    end
end
opts.SelectedVariableNames = intersect(opts.VariableNames, needCols);
Tfull = readtable(csvFile, opts);

% Clamp window to file size
Nfile = height(Tfull);
startId = max(1, startId);
endId   = min(Nfile, endId);
if startId > endId, error('Empty range: startId=%d > endId=%d', startId, endId); end

T = Tfull(startId:endId, :);
n = height(T);
fprintf('Window [%d..%d] -> %d transactions\n', startId, endId, n);
if n==0, error('No transactions in the selected range.'); end

% Pull weight as double (handle string/numeric gracefully)
if ~any(strcmp(T.Properties.VariableNames, wCol)), error('Weight column "%s" not found.', wCol); end
w_raw = T.(wCol);
if iscell(w_raw)
    w = str2double(string(w_raw));
elseif isstring(w_raw)
    w = str2double(w_raw);
elseif isnumeric(w_raw)
    w = double(w_raw);
else
    w = str2double(string(w_raw));
end
w(~isfinite(w)) = 0;     % guard: NaN/Inf -> 0
w = max(0, w(:));        % nonnegative column vector, n-by-1

% Access lists
if ~all(ismember({'access_read','access_write'}, T.Properties.VariableNames))
    error('access_read/access_write columns are required.');
end
Rcol = string(T.access_read);
Wcol = string(T.access_write);

tStart=tic;

%% -------- Parse access lists & build inverted index --------
parse_list = @(s) local_parse_list(s);   % helper below

writers = containers.Map('KeyType','char','ValueType','any');
readers = containers.Map('KeyType','char','ValueType','any');

for i = 1:n
    Rset = parse_list(Rcol(i));
    Wset = parse_list(Wcol(i));

    % index writes
    for t = 1:numel(Wset)
        k = char(Wset(t));
        if isKey(writers,k), writers(k) = [writers(k), i]; else, writers(k) = i; end
    end
    % index reads
    for t = 1:numel(Rset)
        k = char(Rset(t));
        if isKey(readers,k), readers(k) = [readers(k), i]; else, readers(k) = i; end
    end
end

%% -------- Build undirected conflict edges (i<j) --------
% Conflict if they share any key and at least one tx WRITES it: W∩W or W∩R
Ei = []; Ej = [];

wkeys = writers.keys;
for kk = 1:numel(wkeys)
    k = wkeys{kk};
    wlist = writers(k);               % writers of key k

    % writer-writer all pairs
    m = numel(wlist);
    if m >= 2
        for a = 1:m-1
            ia = wlist(a);
            for b = a+1:m
                ib = wlist(b);
                i1 = min(ia,ib); i2 = max(ia,ib);
                Ei(end+1,1) = i1; Ej(end+1,1) = i2; 
            end
        end
    end
    % writer-reader pairs
    if isKey(readers,k)
        rlist = readers(k);
        for a = 1:numel(wlist)
            ia = wlist(a);
            for b = 1:numel(rlist)
                ib = rlist(b);
                if ia==ib, continue; end
                i1 = min(ia,ib); i2 = max(ia,ib);
                Ei(end+1,1) = i1; Ej(end+1,1) = i2; 
            end
        end
    end
end

E = [Ei, Ej];
if ~isempty(E), E = unique(E, 'rows'); end
mEdges = size(E,1);
fprintf('Built conflict graph: n=%d, |E|=%d\n', n, mEdges);

%% -------- Adjacency lists & remaining degrees --------
neighbors = cell(n,1);
if ~isempty(E)
    for ek = 1:size(E,1)
        a = E(ek,1); b = E(ek,2);
        neighbors{a}(end+1) = b;
        neighbors{b}(end+1) = a; 
    end
end
deg = zeros(n,1);
for i = 1:n, deg(i) = numel(neighbors{i}); end

%% -------- Greedy packing by score & key (R rounds, cap p) --------
removed = false(n,1);             % scheduled or not
assign  = zeros(n,1);             % round assignment (0=not selected)
roundSets = cell(R,1);            % selected ids per round
selectedCount = 0;

for roundIdx = 1:R
    % Remaining nodes & order by K(i)=(-s,-w,+deg,+id)
    remIdx = find(~removed).';
    if isempty(remIdx), break; end
    s = w ./ (deg + 1);          % compute score on current graph
    K = [ -s(remIdx), -w(remIdx), deg(remIdx), remIdx(:) ];
    [~, ord] = sortrows(K, [1 2 3 4]);     % deterministic
    order = remIdx(ord);

    % Build round: maximal independent set up to capacity p
    blocked = false(n,1);
    picks = []; picksCap = 0;

    for tIdx = 1:numel(order)
        if picksCap >= p, break; end
        iLoc = order(tIdx);
        if removed(iLoc) || blocked(iLoc), continue; end
        picks(end+1) = iLoc; 
        picksCap = picksCap + 1;
        Ni = neighbors{iLoc};
        if ~isempty(Ni), blocked(Ni) = true; end
    end

    if isempty(picks), break; end

    % Commit picks: remove from graph, update neighbors' degrees
    roundSets{roundIdx} = picks(:);
    assign(picks) = roundIdx;
    removed(picks) = true;
    selectedCount = selectedCount + numel(picks);

    for u = picks
        Nu = neighbors{u};
        for v = Nu
            if ~removed(v)
                deg(v) = deg(v) - 1;
            end
        end
    end
end

usedRounds = find(~cellfun(@isempty, roundSets));
R_used = numel(usedRounds);
fprintf('Selected %d tx across %d/%d rounds (capacity p=%d)\n', selectedCount, R_used, R, p);

%% -------- Build schedule matrix & table (column-wise, robust) --------
xMat = zeros(n, max(R_used,1));
for kRound = 1:numel(usedRounds)
    rId = usedRounds(kRound);
    ids = roundSets{rId};
    xMat(ids, rId) = 1;
end

origId = (startId:endId).';   % map local index -> file row index

% Preallocate columns
Round_col   = zeros(selectedCount,1);
Slot_col    = zeros(selectedCount,1);
LocalTx_col = zeros(selectedCount,1);
FileRow_col = zeros(selectedCount,1);
Weight_col  = zeros(selectedCount,1);

tFill = 0;
for kRound = 1:numel(usedRounds)
    rId = usedRounds(kRound);       % scalar round number
    ids = roundSets{rId};           % vector of local tx indices
    for sidx = 1:numel(ids)
        iLoc = ids(sidx);           % scalar local index in [1..n]
        tFill = tFill + 1;

        Round_col(tFill)   = double(rId);
        Slot_col(tFill)    = double(sidx);
        LocalTx_col(tFill) = double(iLoc);
        FileRow_col(tFill) = double(origId(iLoc));

        ww = w(iLoc);
        if ~isscalar(ww) || ~isfinite(ww), ww = 0; end
        Weight_col(tFill)  = double(ww);
    end
end
if tFill < selectedCount
    Round_col   = Round_col(1:tFill);
    Slot_col    = Slot_col(1:tFill);
    LocalTx_col = LocalTx_col(1:tFill);
    FileRow_col = FileRow_col(1:tFill);
    Weight_col  = Weight_col(1:tFill);
end

Sched = table(Round_col, Slot_col, LocalTx_col, FileRow_col, Weight_col, ...
    'VariableNames', {'Round','Slot','LocalTx','FileRow','Weight'});

disp(Sched(1:min(height(Sched), 40), :));  % show first 40
totalWeight = sum(w(assign>0));
fprintf('Total weight selected = %.6g (sum over assigned tx)\n', totalWeight);
totalSec = toc(tStart);
fprintf('Total time: %.3f s\n', totalSec);

%% -------- Verification --------
tol = 1e-9;
% 1) At-most-one per tx
assignPerTx = sum(xMat,2);
viol_tx = find(assignPerTx > 1 + tol);
% 2) Capacity per round
loadPerRound = sum(xMat,1);
viol_cap = find(loadPerRound > p + tol);
% 3) No conflicts within any round
badPairs = [];   % [round, i, j]
if ~isempty(E)
    for kRound = 1:numel(usedRounds)
        rId = usedRounds(kRound);
        idx = find(xMat(:,rId) > 0.5);
        if numel(idx) > 1
            inRound = false(n,1); inRound(idx) = true;
            hit = inRound(E(:,1)) & inRound(E(:,2));
            if any(hit)
                eidx = find(hit);
                badPairs = [badPairs; [rId*ones(numel(eidx),1), E(eidx,1), E(eidx,2)]]; %#ok<AGROW>
            end
        end
    end
end
% 4) Slots bound
slotsOk = (selectedCount <= R*p + tol);

% Report
if isempty(viol_tx) && isempty(viol_cap) && isempty(badPairs) && slotsOk
    fprintf('Verification: OK  | ≤1 per tx, per-round capacity OK, no conflicts inside rounds.\n');
else
    fprintf('Verification: issues detected \n');
    if ~isempty(viol_tx)
        fprintf('  - At-most-one violated (count=%d). Examples: %s\n', numel(viol_tx), mat2str(viol_tx(1:min(10,end))'));
    end
    if ~isempty(viol_cap)
        fprintf('  - Capacity violated in rounds: %s\n', mat2str(viol_cap));
    end
    if ~isempty(badPairs)
        k = min(10, size(badPairs,1));
        fprintf('  - Conflicts within rounds (showing %d/%d):\n', k, size(badPairs,1));
        disp(array2table(badPairs(1:k,:), 'VariableNames', {'Round','i','j'}));
    end
    if ~slotsOk
        fprintf('  - Selected %d > R*p=%d\n', selectedCount, R*p);
    end
end

%% -------- Local: parse list of keys --------
function toks = local_parse_list(s)
    % Turn forms like:
    %   ["0xA","0xB"]   0xA,0xB   0xA; 0xB   0xA 0xB
    % into a string array of lowercase tokens (may be empty).
    if ismissing(s) || strlength(s)==0
        toks = strings(0,1); return;
    end
    t = lower(strtrim(s));
    % Drop common wrappers and quotes
    t = replace(t, ["[","]","{","}","'",char(34)], "");
    % Normalize separators to space
    t = regexprep(t, '[;,\|\s]+', ' ');
    if strlength(t)==0
        toks = strings(0,1); return;
    end
    toks = string(strsplit(t, ' '));
    toks = toks(toks~="");
end
