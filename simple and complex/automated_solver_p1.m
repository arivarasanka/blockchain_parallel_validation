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

% automated_solver_p1.m
% Batch runner for P1 on a gas-limited slice starting at startId.
% Param CSV must have columns: startId, gasLimit, p
% Each row is executed independently; logs are written per-run.
%
% Requires: solve_min_rounds_p1.m on MATLAB path.

clear; clc;

%% ===== Config =====
paramsCsv = 'p1_batch_params.csv';                 % CSV with startId,gasLimit,p
csvFile   = 'tx_access_21631019_21633879.csv';     % dataset with tx_hash, gas_used, access_read, access_write
logsDir   = 'logs_p1_batch';                       % folder to store per-run logs
maxShow   = 80;                                    % how many scheduled rows to display in each log

% Ensure logs folder exists
if ~exist(logsDir,'dir'), mkdir(logsDir); end

%% ===== Load & validate parameter CSV =====
if ~isfile(paramsCsv), error('Params CSV not found: %s', paramsCsv); end
P = readtable(paramsCsv);

needCols = {'startId','gasLimit','p'};
missing = setdiff(needCols, P.Properties.VariableNames);
if ~isempty(missing)
    error('Params CSV missing required column(s): %s', strjoin(missing, ', '));
end

% Normalize numeric
P.startId  = double(P.startId);
P.gasLimit = double(P.gasLimit);
P.p        = double(P.p);

Njobs = height(P);
fprintf('Loaded %d runs from %s\n', Njobs, paramsCsv);

%% ===== Batch loop with per-run logging & error handling =====
for r = 1:Njobs
    startId  = P.startId(r);
    gasLimit = P.gasLimit(r);
    p        = P.p(r);

    % Unique log filename
    ts  = datestr(now,'yyyymmdd_HHMMSSFFF');
    tag = sprintf('start%d_gas%.0f_p%d', startId, gasLimit, p);
    tag = regexprep(tag, '[^\w\-\.]', '_');
    logf = fullfile(logsDir, sprintf('%s_%s.log', tag, ts));

    fprintf('\n=== RUN %d/%d ===  (startId=%d, gasLimit=%.0f, p=%d)\n', r, Njobs, startId, gasLimit, p);
    fprintf('Log: %s\n', logf);

    diary off;             % reset diary if needed
    diary(logf); diary on; % start per-run logging
    try
        run_one_p1(csvFile, startId, gasLimit, p, maxShow);
        fprintf('RUN %d/%d: SUCCESS\n', r, Njobs);
    catch ME
        fprintf('RUN %d/%d: ERROR\n', r, Njobs);
        fprintf('--- ERROR REPORT START ---\n');
        fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fprintf('--- ERROR REPORT END ---\n');
    end
    diary off; % close this run's log
end

fprintf('\nAll batch runs finished. Logs in: %s\n', logsDir);

%% ================= Helper (single-run) =================
function run_one_p1(csvFile, startId, gasLimit, p, maxShow)

    %% ---------- Read needed columns ----------
    opts = detectImportOptions(csvFile, 'TextType','string');
    needCols = ["tx_hash","gas_used","access_read","access_write"];
    opts.SelectedVariableNames = intersect(opts.VariableNames, needCols);

    % Force text columns to string (robust parsing)
    for c = ["tx_hash","access_read","access_write"]
        if any(strcmp(opts.VariableNames, c))
            idx = find(strcmp(opts.VariableNames, c), 1);
            opts.VariableTypes{idx} = 'string';
        end
    end

    Tfull = readtable(csvFile, opts);
    Ncsv  = height(Tfull);
    if Ncsv == 0, error('No rows in %s', csvFile); end

    % Clamp startId
    startId = max(1, min(startId, Ncsv));

    %% ---------- Build gas-limited window [startId : endId] ----------
    gas_used_col = Tfull.gas_used;
    if isstring(gas_used_col)
        gu_all = str2double(gas_used_col);
    elseif isnumeric(gas_used_col)
        gu_all = double(gas_used_col);
    else
        gu_all = str2double(string(gas_used_col));
    end
    gu_all(~isfinite(gu_all)) = 21000;  % fallback
    gu_all = max(1, round(gu_all));

    gu_tail = gu_all(startId:Ncsv);
    cs = cumsum(gu_tail);
    k = find(cs <= gasLimit, 1, 'last');
    if isempty(k) || k == 0
        error('Row %d alone exceeds gasLimit=%.0f (gas_used=%.0f).', startId, gasLimit, gu_tail(1));
    end
    endId = startId + k - 1;

    % Slice
    T = Tfull(startId:endId, :);
    n = height(T);
    sumGas = sum(gu_all(startId:endId));
    fprintf('Selected rows [%d : %d] of %d. n=%d | sum(gas_used)=%.0f (limit=%.0f)\n', ...
            startId, endId, Ncsv, n, sumGas, gasLimit);
    if n == 0, error('Empty selection after gas slicing.'); end

    %% ---------- Columns & durations ----------
    txhash       = string(T.tx_hash);
    access_read  = string(T.access_read);
    access_write = string(T.access_write);

    if isstring(T.gas_used), t = str2double(T.gas_used); else, t = double(T.gas_used); end
    t(~isfinite(t)) = 21000;
    t = max(1, round(t(:)));

    %% ---------- Parse access lists ----------
    parse_list = @(s) local_parse_list(s);

    writers = containers.Map('KeyType','char','ValueType','any');
    readers = containers.Map('KeyType','char','ValueType','any');

    for i = 1:n
        Rset = parse_list(access_read(i));
        Wset = parse_list(access_write(i));

        for u = 1:numel(Wset)
            k = char(Wset(u));
            if isKey(writers,k), writers(k) = [writers(k), i];
            else,                writers(k) = i;
            end
        end

        for u = 1:numel(Rset)
            k = char(Rset(u));
            if isKey(readers,k), readers(k) = [readers(k), i];
            else,                readers(k) = i;
            end
        end
    end

    %% ---------- Build conflicts (undirected), then direct by i<j ----------
    Ei = []; Ej = [];

    wkeys = writers.keys;
    for kk = 1:numel(wkeys)
        k = wkeys{kk};
        wlist = writers(k);

        % writer-writer
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

        % writer-reader
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

    % For P1: edges already (i<j) by construction
    E = E_und;

    %% ---------- Solve P1 (continuous-time) ----------
    % [xMat, s, e, M, objval, exitflag, output] = solve_min_rounds_p1(...)
    [xMat, s, e, M, objval, exitflag, ~] = solve_min_rounds_p1(n, p, t, E, 'Display','iter');

    fprintf('Makespan M = %.0f\n', M);
    fprintf('exitflag = %d (1 means optimal)\n', exitflag);

    %% ---------- Post-process: core assignment & order ----------
    assignCore = zeros(n,1);
    for i = 1:n
        ci = find(xMat(i,:)>0.5, 1, 'first'); if isempty(ci), ci = 0; end
        assignCore(i) = ci;
    end

    roundInCore = zeros(n,1);
    for c = 1:p
        coreJobs = find(xMat(:,c) > 0.5);
        if isempty(coreJobs), continue; end
        [~, ksort] = sortrows([s(coreJobs), e(coreJobs), coreJobs], [1 1 1]);
        coreJobs = coreJobs(ksort);
        roundInCore(coreJobs) = 1:numel(coreJobs);
    end

    %% ---------- Print compact schedule ----------
    csv_rows = (startId:(startId+n-1)).';
    S = table(csv_rows, assignCore, roundInCore, t, ...
        'VariableNames', {'csv_row','Core','Order','t'});
    if any(txhash ~= "")
        S = [table(txhash,'VariableNames',{'tx_hash'}) S];
    end

    Nprint = min(maxShow, height(S));
    fprintf('\n--- Schedule (first %d rows) ---\n', Nprint);
    disp(S(1:Nprint,:));

    %% ---------- Local parser (nested) ----------
    function toks = local_parse_list(s)
        % Robustly parse access list cell into a vector of lowercased tokens.
        if ismissing(s) || strlength(s) == 0
            toks = strings(0,1); return;
        end
        tloc = lower(strtrim(s));
        tloc = replace(tloc, ["[","]","{","}","'",char(34)], "");
        tloc = regexprep(tloc, '[;,\|\s]+', ' ');
        if strlength(tloc) == 0
            toks = strings(0,1); return;
        end
        toks = string(strsplit(tloc, ' '));
        toks = toks(toks ~= "");
    end
end
