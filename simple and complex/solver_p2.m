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

% solver_p2.m — P2 (complex) with conflicts from access lists
% Uses:
%   t_i  <- gas_used
%   w_i  <- weight_theoretical_eth
% Conflicts derived from access_read / access_write (no E.csv).
% NEW: select a contiguous CSV row range [rowFrom : rowTo] before solving.

clear; clc;

% ---- Config ----
p        = 8;                                        % #cores
R        = 425280;                                       % horizon; [] => sum(t)
txsCsv   = 'tx_access_21631019_21631879.csv';        % headered CSV
rowFrom  = 1;                                        % first CSV row to include (1-based)
rowTo    = 60;                                      % last CSV row to include (inclusive)
% -----------------------

% ---- Read CSV (only needed columns) ----
opts = detectImportOptions(txsCsv, 'TextType','string');
need = ["gas_used","weight_theoretical_wei","access_read","access_write"];
opts.SelectedVariableNames = intersect(opts.VariableNames, need);

% Force text columns to string so parsing is simple
for c = ["access_read","access_write"]
    if any(strcmp(opts.VariableNames, c))
        idx = find(strcmp(opts.VariableNames, c), 1);
        opts.VariableTypes{idx} = 'string';
    end
end

Tall = readtable(txsCsv, opts);
Ncsv = height(Tall);
if Ncsv==0, error('No rows in %s', txsCsv); end

% Clamp & subset
rowFrom = max(1, min(rowFrom, Ncsv));
rowTo   = max(rowFrom, min(rowTo, Ncsv));
fprintf('Selected CSV rows [%d : %d] of %d total.\n', rowFrom, rowTo, Ncsv);

T = Tall(rowFrom:rowTo, :);
n = height(T);

% t <- gas_used
if any(strcmp(T.Properties.VariableNames,'gas_used'))
    if isstring(T.gas_used), t = str2double(T.gas_used); else, t = double(T.gas_used); end
else
    error('Column "gas_used" not found in %s', txsCsv);
end
t(~isfinite(t)) = 21000;                    % fallback
t = max(1, round(t(:)));

% w <- weight_theoretical_eth
if any(strcmp(T.Properties.VariableNames,'weight_theoretical_wei'))
    if isstring(T.weight_theoretical_wei), w = str2double(T.weight_theoretical_wei);
    else, w = double(T.weight_theoretical_wei); end
else
    error('Column "weight_theoretical_wei" not found in %s', txsCsv);
end
w(~isfinite(w)) = 0;
w = w(:);

fprintf('n(subset)=%d | using gas_used as t and weight_theoretical_wei as w\n', n);
%%

% ---- Build conflicts from access_read / access_write ----
access_read  = strings(n,1);
access_write = strings(n,1);
if any(strcmp(T.Properties.VariableNames,'access_read')),  access_read  = string(T.access_read);  end
if any(strcmp(T.Properties.VariableNames,'access_write')), access_write = string(T.access_write); end
access_read(ismissing(access_read))   = "";
access_write(ismissing(access_write)) = "";

% Parse lists and index readers/writers by key
writers = containers.Map('KeyType','char','ValueType','any');
readers = containers.Map('KeyType','char','ValueType','any');

parse_list = @(s) local_parse_list(s);

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

% Undirected conflicts: (W∩W) ∪ (W∩R) ∪ (R∩W). (Pure R∩R is NOT a conflict.)
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

E = [Ei, Ej];
if ~isempty(E), E = unique(E, 'rows'); end
fprintf('Conflicts built from access lists: |E|=%d\n', size(E,1));

% ---- Solve ----
if isempty(R)
    [xMat, s, e, objval, exitflag, output] = solve_min_rounds_p2(n, p, t, E, w, 'Display','iter');
else
    [xMat, s, e, objval, exitflag, output] = solve_min_rounds_p2(n, p, t, E, w, 'R', R, 'Display','iter');
end

% ---- Postprocess ----
selected = any(xMat > 0.5, 2);
core     = zeros(n,1);
for i = 1:n
    if selected(i), [~, core(i)] = max(xMat(i,:)); end
end
M = any(selected) * max(e(selected));

Sched = table((rowFrom:rowTo).', (1:n).', core, s, e, t(:), w(:), selected, ...
    'VariableNames', {'csv_row','local_id','core','s','e','t','w','selected'});

% ---- Summary ----
fprintf('\n=== P3 (complex) Solve Summary ===\n');
fprintf('rows [%d:%d] → n=%d, p=%d, selected=%d, objective=%.6g, makespan=%.0f\n', ...
    rowFrom, rowTo, n, p, sum(selected), objval, M);
fprintf('exitflag=%d\n', exitflag);

% Show first 20 selected rows
selIdx = find(selected);
k = min(100, numel(selIdx));
if k > 0
    fprintf('\nSelected (first %d rows):\n', k);
    disp(Sched(selIdx(1:k), :));
else
    fprintf('\nNo transactions selected.\n');
end

% ===== Local parser =====
function toks = local_parse_list(s)
    % Parse an access list cell into a vector of lowercased tokens.
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
