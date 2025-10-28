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

%% heuristic_p1_complex.m
% Problem 1 heuristic (height + volume + outdegree, tie by id) on a CSV row range.
% Execution time of each transaction i is t_i = gas_used(i).

clear; clc;

%% ---------- Config ----------
csvFile = 'tx_access_21631019_21633879.csv';  % CSV with: tx_hash, access_read, access_write, gas_used, (optional) transactionIndex
startId =4290;                                   % first CSV row to include (1-based)
endId   =4299;                                 % last CSV row to include (inclusive)
p       = 2;                                   % number of processors (parallel slots)
% --------------------------------

%% ---------- Read needed columns & slice by row range ----------
opts = detectImportOptions(csvFile, 'TextType','string');

needCols = ["tx_hash","access_read","access_write","gas_used","transactionIndex"];
missing  = needCols(~ismember(needCols, string(opts.VariableNames)));
% transactionIndex is optional; others must exist
if ~isempty(setdiff(missing,"transactionIndex"))
    error('Missing required column(s): %s', strjoin(cellstr(setdiff(missing,"transactionIndex")), ', '));
end

% Force string for text columns
for c = ["tx_hash","access_read","access_write"]
    idx = find(strcmp(opts.VariableNames, c), 1);
    if ~isempty(idx), opts.VariableTypes{idx} = 'string'; end
end
% Ensure numerics
for c = ["gas_used","transactionIndex"]
    if any(strcmp(opts.VariableNames, c))
        idx = find(strcmp(opts.VariableNames, c), 1);
        if ~isempty(idx), opts.VariableTypes{idx} = 'double'; end
    end
end

opts.SelectedVariableNames = intersect(cellstr(needCols), opts.VariableNames);
Tfull = readtable(csvFile, opts);
Ncsv  = height(Tfull);
if Ncsv == 0, error('No rows in %s', csvFile); end

% Clamp range
startId = max(1, min(startId, Ncsv));
endId   = max(1, min(endId,   Ncsv));
if endId < startId, tmp = startId; startId = endId; endId = tmp; end

% Slice
T = Tfull(startId:endId, :);
n = height(T);
fprintf('Selected CSV rows [%d : %d] of %d total ? %d transactions.\n', startId, endId, Ncsv, n);
if n == 0, error('Empty selection.'); end

% Columns
tx_hash      = string(T.tx_hash);
access_read  = string(T.access_read);
access_write = string(T.access_write);

% Execution times (weights): t_i = gas_used
if any(strcmp(T.Properties.VariableNames,"gas_used"))
    t = double(T.gas_used);
else
    error('gas_used column is required.');
end
t(~isfinite(t)) = 0;       % sanitize NaN/Inf ? 0
t(t < 0) = 0;              % no negatives
t = t(:);

% Deterministic id(i) for tie-breaking:
% prefer transactionIndex if present; else use local row order.
if any(strcmp(T.Properties.VariableNames,"transactionIndex"))
    id = double(T.transactionIndex);
    if any(~isfinite(id)), id = (1:n)'; end
else
    id = (1:n)';   % local order
end

% Helper: normalize & tokenize access lists
parse_list = @local_parse_list;
tstart=tic;

%% ---------- Build UNDIRECTED conflict graph from R/W sets ----------
% Conflict iff two txs share a key and at least one WRITES that key:
%   WW or WR or RW   (RR alone is NOT a conflict).
writers = containers.Map('KeyType','char','ValueType','any');  % key -> [local idx ...]
readers = containers.Map('KeyType','char','ValueType','any');

for i = 1:n
    Rset = parse_list(access_read(i));
    Wset = parse_list(access_write(i));

    % index writes
    for tt = 1:numel(Wset)
        k = char(Wset(tt));
        if isKey(writers,k), writers(k) = [writers(k), i];
        else,                writers(k) = i;
        end
    end
    % index reads
    for tt = 1:numel(Rset)
        k = char(Rset(tt));
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

%% ---------- Convert conflicts to precedence arcs (acyclic) ----------
% Prefer orientation by transactionIndex (earlier -> later). If missing, E_und already i<j by row.
if any(strcmp(T.Properties.VariableNames,"transactionIndex"))
    % position in order by transactionIndex (stable)
    [~, orderByTxIdx] = sortrows([id, (1:n)'], [1 2]);
    pos = zeros(n,1); pos(orderByTxIdx) = 1:n;
    if ~isempty(E_und)
        iPos = pos(E_und(:,1));
        jPos = pos(E_und(:,2));
        mask_ij = iPos <= jPos;
        E = [ E_und(mask_ij,1), E_und(mask_ij,2); ...
              E_und(~mask_ij,2), E_und(~mask_ij,1) ];
    else
        E = zeros(0,2);
    end
else
    % E_und already built with (i<j) by construction
    E = E_und;
    pos = (1:n)';   % row order as "position"
end
m = size(E,1);
fprintf('Precedence arcs built: %d\n', m);

%% ---------- Build succ/pred/indeg/outdeg in ONE scan ----------
succ   = cell(n,1);
pred   = cell(n,1);
indeg  = zeros(n,1);
outdeg = zeros(n,1);

for k = 1:m
    i = E(k,1); j = E(k,2);
    succ{i}(end+1) = j;        
    pred{j}(end+1) = i;        
    indeg(j)  = indeg(j)  + 1;
    outdeg(i) = outdeg(i) + 1;
end

%% ---------- Heights h(i) & Volumes v(i) (weighted) ----------
% h(i) = t_i + max_{j in succ(i)} h(j); sinks have h = t_i
% v(i) = t_i + sum_{j in succ(i)} v(j); sinks have v = t_i (double-counts merges)
h = t;
v = t;

% Backward sweep in a *topological* descending order.
% Use 'pos' (earlier -> later) so descending pos processes successors first.
[~, topo_desc] = sort(pos, 'descend');
for uu = topo_desc(:).'
    Si = succ{uu};
    if ~isempty(Si)
        h(uu) = t(uu) + max(h(Si));
        v(uu) = t(uu) + sum(v(Si));
    else
        h(uu) = t(uu);
        v(uu) = t(uu);
    end
end

%% ---------- Heuristic scheduling (event-driven, p processors) ----------
% Key K(i) = (-h(i), -v(i), -outdeg(i), +id(i)); pick lexicographically smallest
% when assigning from ready.

% State
remaining = true(n,1);                % not yet finished
indeg_cur = indeg;                    % mutable indegrees
ready     = find(indeg_cur == 0).';   % row vector
running   = false(n,1);               % currently running?
startTime = NaN(n,1);
finishTime= NaN(n,1);
whichCPU  = NaN(n,1);

timeNow   = 0;
freeCPUs  = 1:p;                      % free processor ids
numRunning= 0;

% Initialize
[ready, running, startTime, finishTime, whichCPU, freeCPUs, numRunning] = ...
    assign_from_ready(ready, freeCPUs, h, v, outdeg, id, ...
                      running, startTime, finishTime, whichCPU, timeNow, t, numRunning);

% Main event loop
while any(remaining)
    if numRunning == 0
        % no job is running but some remain -> must have something ready
        if isempty(ready)
            error('Deadlock detected (cycle). Check precedence construction.');
        end
        [ready, running, startTime, finishTime, whichCPU, freeCPUs, numRunning] = ...
            assign_from_ready(ready, freeCPUs, h, v, outdeg, id, ...
                              running, startTime, finishTime, whichCPU, timeNow, t, numRunning);
        if numRunning == 0
            error('Unable to start any job though ready is nonempty.');
        end
    end

    % Advance time to next finish event
    activeIdx   = find(running);
    [tNext, kmin] = min(finishTime(activeIdx));
    jobFin      = activeIdx(kmin);
    timeNow     = tNext;

    % Complete jobFin
    running(jobFin)    = false;
    remaining(jobFin)  = false;
    numRunning         = numRunning - 1;
    % free its CPU
    freeCPUs(end+1) = whichCPU(jobFin); 

    % Update indegrees of successors; enqueue any that become ready
    Sj = succ{jobFin};
    for tt2 = 1:numel(Sj)
        vtx = Sj(tt2);
        indeg_cur(vtx) = indeg_cur(vtx) - 1;
        if indeg_cur(vtx) == 0
            ready(end+1) = vtx; 
        end
    end

    % Fill idle CPUs from ready set
    [ready, running, startTime, finishTime, whichCPU, freeCPUs, numRunning] = ...
        assign_from_ready(ready, freeCPUs, h, v, outdeg, id, ...
                          running, startTime, finishTime, whichCPU, timeNow, t, numRunning);
end

Cmax = max(finishTime);
fprintf('Heuristic finished: makespan = %.0f\n', Cmax);
totalSec = toc(tstart);
fprintf('Total time: %.3f s\n', totalSec);

%% ---------- Summaries (compact) ----------
S = table((startId:endId).', id, startTime, finishTime, whichCPU, ...
    'VariableNames', {'csv_row','id','start','finish','cpu'});
S = sortrows(S, {'cpu','start','id'});
disp(S)

%% ---------- Verification: constraints & sanity checks ----------
tol = 1e-9; ok = true;

% 1) Every job scheduled with finite start/finish & assigned to exactly one CPU
if any(~isfinite(startTime)) || any(~isfinite(finishTime))
    ok = false; fprintf('? Some jobs have NaN start/finish.\n');
end
if any(isnan(whichCPU))
    ok = false; fprintf('? Some jobs have no CPU assigned.\n');
end
xMat = zeros(n, p);
for i = 1:n
    if isfinite(whichCPU(i)) && whichCPU(i)>=1
        xMat(i, whichCPU(i)) = 1;
    end
end
rowSums = sum(xMat,2);
bad = find(abs(rowSums - 1) > tol);
if ~isempty(bad)
    ok = false;
    fprintf('? Exactly-one violated for %d jobs. Example: %s\n', ...
        numel(bad), mat2str(bad(1:min(10,numel(bad)))'));
end

% 2) Duration consistency: e - s == t
dur = finishTime - startTime;
bad = find(abs(dur - t) > tol);
if ~isempty(bad)
    ok = false;
    fprintf('? Duration mismatch for %d jobs. Example: %s\n', ...
        numel(bad), mat2str(bad(1:min(10,numel(bad)))'));
end

% 3) Precedence: start(j) >= finish(i) for all (i -> j)
prec_bad = [];
if ~isempty(E)
    si = startTime(E(:,2));
    ei = finishTime(E(:,1));
    bad_mask = (si + tol < ei);
    if any(bad_mask)
        idx = find(bad_mask);
        prec_bad = [E(idx,1:2), ei(idx), si(idx)];  % [i, j, e(i), s(j)]
        ok = false;
        k = min(10, size(prec_bad,1));
        fprintf('? Precedence violations (showing %d/%d): [i j e(i) s(j)]\n', k, size(prec_bad,1));
        disp(array2table(prec_bad(1:k,:), 'VariableNames', {'i','j','e_i','s_j'}));
    end
end

% 4) No overlap on any core
overlap_cnt = 0; sample_ov = [];
for c = 1:p
    jobs = find(xMat(:,c) > 0.5);
    if numel(jobs) <= 1, continue; end
    [~,ord] = sort(startTime(jobs)); jobs = jobs(ord);
    for k = 1:numel(jobs)-1
        i1 = jobs(k); i2 = jobs(k+1);
        if startTime(i2) < finishTime(i1) - tol
            overlap_cnt = overlap_cnt + 1;
            if size(sample_ov,1) < 10
                sample_ov(end+1,:) = [c, i1, i2, finishTime(i1), startTime(i2)]; 
            end
        end
    end
end
if overlap_cnt > 0
    ok = false;
    fprintf('? Overlaps on cores: %d cases. Examples [core i j e(i) s(j)]:\n', overlap_cnt);
    disp(array2table(sample_ov, 'VariableNames', {'core','i','j','e_i','s_j'}));
end

% 5) Global capacity check (<= p running at any time)
events = [startTime,  ones(n,1); ...
          finishTime, -ones(n,1)];
events = sortrows(events, [1 2]);   % finishes (-1) before starts (+1) at same time
active = 0; cap_violate = false; maxActive = 0;
for u = 1:size(events,1)
    active = active + events(u,2);
    if active > maxActive, maxActive = active; end
    if active > p + tol
        cap_violate = true; break;
    end
end
if cap_violate
    ok = false;
    fprintf('? Capacity violation: concurrent=%d > p=%d\n', active, p);
else
    fprintf('Max concurrent jobs observed: %d (p=%d)\n', maxActive, p);
end

% 6) Makespan consistency
Mcalc = max(finishTime);
if abs(Cmax - Mcalc) > tol
    ok = false;
    fprintf('? Makespan mismatch: Cmax=%.9g vs max(finishTime)=%.9g\n', Cmax, Mcalc);
end

if ok
    fprintf('? Verification: OK ? exactly-one, duration, precedence, no-overlap, capacity, makespan.\n');
end

%% ---------- Local functions ----------
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

function ord = order_ready(ready, h, v, outdeg, id)
    % Return ordering of 'ready' by K(i) = (-h, -v, -outdeg, +id)
    if isempty(ready), ord = []; return; end
    Kr = [ -h(ready), -v(ready), -outdeg(ready), id(ready) ];
    [~, ord] = sortrows(Kr, [1 2 3 4]);  % lexicographic ascending
end

function [ready, running, startTime, finishTime, whichCPU, freeCPUs, numRunning] = ...
    assign_from_ready(ready, freeCPUs, h, v, outdeg, id, ...
                      running, startTime, finishTime, whichCPU, timeNow, t, numRunning)
    % Assign as many ready jobs as there are free CPUs, using the K-order.
    if isempty(ready) || isempty(freeCPUs), return; end
    ord = order_ready(ready, h, v, outdeg, id);
    ready = ready(ord);  % smallest K first
    take = min(numel(ready), numel(freeCPUs));
    picked = ready(1:take);
    ready(1:take) = [];
    for uu = picked
        cpu = freeCPUs(1);
        freeCPUs(1) = [];
        whichCPU(uu)  = cpu;
        startTime(uu) = timeNow;
        finishTime(uu)= timeNow + t(uu);
        running(uu)   = true;
        numRunning    = numRunning + 1;
    end
end
