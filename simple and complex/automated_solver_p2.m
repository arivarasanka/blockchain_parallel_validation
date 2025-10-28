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

% automated_solver_p2.m
% Batch runner for P2 (complex) using conflicts from access lists.
% Reads a parameter CSV with columns: p, R, x, rowFrom
% For each row, runs the solver and writes a separate log file.
% Errors are caught and logged; the batch continues.

clear; clc;

% ==== Config ====
paramsCsv = 'batch_params.csv';                 % CSV with columns: p,R,x,rowFrom
txsCsv    = 'tx_access_21631019_21633879.csv';  % The big dataset to read transactions from
logsDir   = 'logs_p3_batch';                    % Folder to store per-run logs
maxShowSelected = 100;                          % Max rows to print from selected schedule

% Ensure logs folder exists
if ~exist(logsDir,'dir'), mkdir(logsDir); end

% ---- Load and validate parameter CSV ----
if ~isfile(paramsCsv), error('Params CSV not found: %s', paramsCsv); end
P = readtable(paramsCsv);

needCols = {'p','R','x','rowFrom'};
missing = setdiff(needCols, P.Properties.VariableNames);
if ~isempty(missing)
    error('Params CSV missing required column(s): %s', strjoin(missing, ', '));
end

% Normalize/validate types
P.p       = double(P.p);
P.R       = double(P.R);
P.x       = double(P.x);
P.rowFrom = double(P.rowFrom);

Njobs = height(P);
fprintf('Loaded %d runs from %s\n', Njobs, paramsCsv);

% ---- Batch loop with per-run logging and error handling ----
for r = 1:Njobs
    p       = P.p(r);
    R       = P.R(r);
    x       = P.x(r);
    rowFrom = P.rowFrom(r);

    % Build a stable, unique log filename
    ts   = datestr(now,'yyyymmdd_HHMMSSFFF');
    tag  = sprintf('p%d_R%g_x%g_from%d', p, R, x, rowFrom);
    tag  = regexprep(tag, '[^\w\-\.]', '_');  % sanitize
    logf = fullfile(logsDir, sprintf('%s_%s.log', tag, ts));

    fprintf('\n=== RUN %d/%d ===\n', r, Njobs);
    fprintf('Log: %s\n', logf);

    % Start logging
    diary off;  % in case a previous run crashed with diary on
    diary(logf);
    diary on;

    try
        run_one(p, R, x, rowFrom, txsCsv, maxShowSelected);
        fprintf('RUN %d/%d: SUCCESS\n', r, Njobs);
    catch ME
        fprintf('RUN %d/%d: ERROR\n', r, Njobs);
        fprintf('--- ERROR REPORT START ---\n');
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fprintf('--- ERROR REPORT END ---\n');
    end

    % Stop logging for this run
    diary off;
end

fprintf('\nAll batch runs finished. Logs in: %s\n', logsDir);

% ====================== Helper: run_one ======================
function run_one(p, R, x, rowFrom, txsCsv, maxShowSelected)
    % One full execution (your single-run logic), slightly modularized.

    fprintf('Params: p=%d, R=%g, x=%g, rowFrom=%d\n', p, R, x, rowFrom);
    fprintf('Source CSV: %s\n', txsCsv);

    % ---- Read CSV (only needed columns) ----
    opts = detectImportOptions(txsCsv, 'TextType','string');
    need = ["gas_used","weight_theoretical_eth","access_read","access_write"];
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

    % ---- Prepare gas_used as numeric for cumulative stopping ----
    if ~any(strcmp(Tall.Properties.VariableNames,'gas_used'))
        error('Column "gas_used" not found in %s', txsCsv);
    end
    if isstring(Tall.gas_used)
        gas_all = str2double(Tall.gas_used);
    else
        gas_all = double(Tall.gas_used);
    end
    gas_all(~isfinite(gas_all)) = 21000;          % fallback for bad/missing cells
    gas_all = max(1, round(gas_all(:)));          % ensure positive integer-ish

    % ---- Decide rowTo so that sum(gas_used) <= x*R ----
    rowFrom = max(1, min(rowFrom, Ncsv));
    target  = p * x * R;
    
    cs = cumsum(gas_all(rowFrom:end));
    
    % choose the largest k with cs(k) <= target
    k = find(cs <= target, 1, 'last');
    
    if isempty(k)
        error('Row %d alone exceeds target=%.0f (gas_used=%g).', ...
              rowFrom, target, gas_all(rowFrom));
    end
    
    rowTo   = rowFrom + k - 1;
    sumGas  = cs(k);
    fprintf('Selected rows [%d:%d]; sum(gas_used)=%g (cap=%g = x*R)\n', ...
            rowFrom, rowTo, sumGas, target);


    % ---- Subset table ----
    T = Tall(rowFrom:rowTo, :);
    n = height(T);

    % t <- gas_used (for the subset)
    if isstring(T.gas_used), t = str2double(T.gas_used); else, t = double(T.gas_used); end
    t(~isfinite(t)) = 21000;
    t = max(1, round(t(:)));

    % w <- weight_theoretical_eth
    if any(strcmp(T.Properties.VariableNames,'weight_theoretical_eth'))
        if isstring(T.weight_theoretical_eth), w = str2double(T.weight_theoretical_eth);
        else, w = double(T.weight_theoretical_eth); end
    else
        error('Column "weight_theoretical_eth" not found in %s', txsCsv);
    end
    w(~isfinite(w)) = 0;
    w = w(:);

    fprintf('n(selected)=%d | using gas_used as t and weight_theoretical_eth as w\n', n);

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

    % Undirected conflicts: (W?W) ? (W?R) ? (R?W). (Pure R?R is NOT a conflict.)
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

    E = [Ei, Ej];
    if ~isempty(E), E = unique(E, 'rows'); end
    fprintf('Conflicts built from access lists: |E|=%d\n', size(E,1));

    % ---- Solve ----
    % NOTE: solve_min_rounds_p2 must be on path.
    [xMat, s, e, objval, exitflag, output] = solve_min_rounds_p2(n, p, t, E, w, 'R', R, 'Display','iter');

    % ---- Postprocess ----
    selected = any(xMat > 0.5, 2);
    core     = zeros(n,1);
    for i = 1:n
        if selected(i), [~, core(i)] = max(xMat(i,:)); end
    end
    M = any(selected) * max(e(selected));

    Sched = table((rowFrom:(rowFrom+n-1)).', (1:n).', core, s, e, t(:), w(:), selected, ...
        'VariableNames', {'csv_row','local_id','core','s','e','t','w','selected'});

    % ---- Summary ----
    fprintf('\n=== P2 (complex) Solve Summary ===\n');
    fprintf('rows [%d:%d] ? n=%d, p=%d, selected=%d, objective=%.6g, makespan=%.0f\n', ...
        rowFrom, rowFrom+n-1, n, p, sum(selected), objval, M);
    fprintf('exitflag=%d\n', exitflag);

    % Show first up to maxShowSelected selected rows
    selIdx = find(selected);
    k = min(maxShowSelected, numel(selIdx));
    if k > 0
        fprintf('\nSelected (first %d rows):\n', k);
        disp(Sched(selIdx(1:k), :));
    else
        fprintf('\nNo transactions selected.\n');
    end
end

% ===== Local parser =====
function toks = local_parse_list(s)
    % Parse an access list cell into a vector of lowercased tokens.
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
