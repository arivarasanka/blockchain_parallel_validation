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

%% solver_p1.m : P1 on a gas-limited slice starting at startId
clear; clc;

%% ---------- Config ----------
csvFile = 'tx_access_21631019_21633879.csv';  % has headers incl. access_read, access_write
startId = 4290;                                   % first CSV row to include (1-based)
gasLimit = 1260000;                                  % stop BEFORE exceeding this total gas
p       = 2;                                    % number of cores

%% ---------- Read needed columns (simple & robust) ----------
opts = detectImportOptions(csvFile, 'TextType','string');

needCols = ["tx_hash","gas_used","access_read","access_write"];
for c = needCols
    if any(strcmp(opts.VariableNames, c))
        if c ~= "gas_used"
            idx = find(strcmp(opts.VariableNames, c), 1);
            opts.VariableTypes{idx} = 'string';   % text as string
        end
    end
end
opts.SelectedVariableNames = intersect(opts.VariableNames, needCols);

Tfull = readtable(csvFile, opts);
Ncsv  = height(Tfull);
if Ncsv == 0, error('No rows in %s', csvFile); end

% Clamp startId
startId = max(1, min(startId, Ncsv));

%% ---------- Build gas-limited contiguous window [startId : endId] ----------
% Pull gas_used robustly from startId..Ncsv
gas_used_col = Tfull.gas_used;
if isstring(gas_used_col)
    gu_all = str2double(gas_used_col);
elseif isnumeric(gas_used_col)
    gu_all = double(gas_used_col);
else
    gu_all = str2double(string(gas_used_col));
end
gu_all(~isfinite(gu_all)) = 21000;         % fallback for bad/empty
gu_all = max(1, round(gu_all));

% Work from startId to the end, accumulate until about to exceed gasLimit
gu_tail = gu_all(startId:Ncsv);
cs = cumsum(gu_tail);
k = find(cs <= gasLimit, 1, 'last');       % number of rows we can take without exceeding
if isempty(k) || k == 0
    error('First transaction at row %d already exceeds gasLimit=%.0f (gas_used=%.0f).', ...
          startId, gasLimit, gu_tail(1));
end
endId = startId + k - 1;

% Slice the chosen window
T = Tfull(startId:endId, :);
n = height(T);
sumGas = sum(gu_all(startId:endId));
fprintf('Selected CSV rows [%d : %d] of %d total. Transactions: %d | sum(gas_used)=%.0f (limit=%.0f)\n', ...
        startId, endId, Ncsv, n, sumGas, gasLimit);
if n == 0, error('Empty selection.'); end

%% ---------- Pull columns ----------
txhash       = string(T.tx_hash);
access_read  = string(T.access_read);
access_write = string(T.access_write);

% Durations from gas_used (already validated)
if isstring(T.gas_used), t = str2double(T.gas_used); else, t = double(T.gas_used); end
t(~isfinite(t)) = 21000;
t = max(1, round(t(:)));

%% ---------- Parse access lists into tokens ----------
parse_list = @(s) local_parse_list(s);

% Build reader/writer indexes by key
writers = containers.Map('KeyType','char','ValueType','any');
readers = containers.Map('KeyType','char','ValueType','any');

for i = 1:n
    Rset = parse_list(access_read(i));
    Wset = parse_list(access_write(i));

    for tkn = 1:numel(Wset)
        k = char(Wset(tkn));
        if isKey(writers,k), writers(k) = [writers(k), i];
        else,                writers(k) = i;
        end
    end

    for tkn = 1:numel(Rset)
        k = char(Rset(tkn));
        if isKey(readers,k), readers(k) = [readers(k), i];
        else,                readers(k) = i;
        end
    end
end

%% ---------- Build conflict graph (undirected), then direct by i<j ----------
Ei = []; Ej = [];

wkeys = writers.keys;
for kk = 1:numel(wkeys)
    k = wkeys{kk};
    wlist = writers(k);

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

E_und = [Ei, Ej];
if ~isempty(E_und), E_und = unique(E_und, 'rows'); end
fprintf('Conflicts (undirected) built: |E| = %d\n', size(E_und,1));

% For P1: edges are already with (i<j) by min/max above
E = E_und;

%% ---------- Solve P1 (continuous-time) ----------
% Assumes solve_min_rounds_p1 is on path and returns:
% [xMat, s, e, M, objval, exitflag, output]
[xMat, s, e, M, objval, exitflag, output] = solve_min_rounds_p1(n, p, t, E, ...
    'Display','iter'); 

fprintf('Makespan M = %.0f\n', M);
fprintf('exitflag = %d (1 means optimal)\n', exitflag);

%% ---------- Post-process: core + within-core order ----------
assignCore = zeros(n,1);
for i = 1:n
    ci = find(xMat(i,:)>0.5, 1, 'first'); if isempty(ci), ci = 0; end
    assignCore(i) = ci;
end

roundInCore = zeros(n,1);
for c = 1:p
    coreJobs = find(xMat(:,c) > 0.5);
    if isempty(coreJobs), continue; end
    % Order by start, then end, then id (deterministic)
    [~, ksort] = sortrows([s(coreJobs), e(coreJobs), coreJobs], [1 1 1]);
    coreJobs = coreJobs(ksort);
    roundInCore(coreJobs) = 1:numel(coreJobs);
end

% Print CSV row, tx_hash (if present), Core, Order, t
csv_rows = (startId:endId).';
if any(txhash ~= "")
    disp(table(csv_rows, txhash, assignCore, roundInCore, t, ...
        'VariableNames', {'csv_row','tx_hash','Core','Order','t'}));
else
    disp(table(csv_rows, assignCore, roundInCore, t, ...
        'VariableNames', {'csv_row','Core','Order','t'}));
end

%% ---------- Local parser ----------
function toks = local_parse_list(s)
    % Robustly parse access list cell into a vector of lowercased tokens.
    if ismissing(s) || strlength(s) == 0
        toks = strings(0,1); return;
    end
    t = lower(strtrim(s));
    t = replace(t, ["[","]","{","}","'",char(34)], "");
    t = regexprep(t, '[;,\|\s]+', ' ');
    if strlength(t) == 0
        toks = strings(0,1); return;
    end
    toks = string(strsplit(t, ' '));
    toks = toks(toks ~= "");
end
