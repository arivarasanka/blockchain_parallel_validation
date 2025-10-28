%% automated_solver_p2.m — Batch P2 on tx-id slices, logging per range
clear; clc;

%% -------- Static config --------
csvFile = 'tx_access_simple_21631019_21635079.csv';
wCol    = 'weight_theoretical_eth';  % weight column to use
rangesFile = 'ranges.csv';           % CSV with headers: startId,endId,p,R
logsDir    = 'logs';                 % where to store per-iteration logs
% -------------------------------------------

%% -------- Load ranges --------
if ~isfile(rangesFile)
    error('Ranges file not found: %s', rangesFile);
end
TR = readtable(rangesFile);
needCols = {'startId','endId','p','R'};
if ~all(ismember(needCols, TR.Properties.VariableNames))
    error('ranges.csv must contain columns: startId,endId,p,R');
end

% Make sure logs folder exists
if ~exist(logsDir, 'dir'), mkdir(logsDir); end

%% -------- Read needed columns once --------
opts = detectImportOptions(csvFile, 'TextType','string');
needColsMain = ["tx_hash","blockNumber","access_read","access_write",wCol];
have = ismember(needColsMain, string(opts.VariableNames));
if ~all(have)
    missing = needColsMain(~have);
    error('Missing required column(s) in CSV: %s', strjoin(cellstr(missing), ', '));
end
for c = ["tx_hash","access_read","access_write"]
    idx = find(strcmp(opts.VariableNames, c), 1);
    opts.VariableTypes{idx} = 'string';
end
opts.SelectedVariableNames = cellstr(needColsMain);
Tfull = readtable(csvFile, opts);
Nfull = height(Tfull);

%% -------- Iterate over ranges and run your solver code --------
for kRun = 1:height(TR)
    % Pull params for this run
    startId = TR.startId(kRun);
    endId   = TR.endId(kRun);
    p       = TR.p(kRun);
    R       = TR.R(kRun);

    % Make a log file name
    stamp = datestr(now,'yyyymmdd_HHMMSS_FFF');
    logName = sprintf('p2_%d_%d_p%d_R%d_%s.log', startId, endId, p, R, stamp);
    logPath = fullfile(logsDir, logName);

    % Start logging this run
    diary off;  % ensure clean
    diary(logPath);
    fprintf('===== Run %d/%d | startId=%d, endId=%d, p=%d, R=%d =====\n', ...
        kRun, height(TR), startId, endId, p, R);

    try
        %% -------- Bound and slice by row ids --------
        sId = max(1, min(startId, Nfull));
        eId = max(1, min(endId,   Nfull));
        if eId < sId
            error('endId (%d) < startId (%d).', eId, sId);
        end
        T = Tfull(sId:eId, :);
        n = height(T);

        fprintf('Slice [%d..%d] → %d transactions.\n', sId, eId, n);
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
        for kk2 = 1:numel(wkeys)
            k = wkeys{kk2};
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

        %% -------- Build in-memory schedule table --------
        rows = R * p;
        R_col = zeros(rows,1);
        P_col = zeros(rows,1);
        IXcol = zeros(rows,1);      % local row index (1..n)
        W_col = zeros(rows,1);
        H_col = strings(rows,1);    % tx_hash (if available)

        t = 0;
        for r = 1:R
            idx = find(xMat(:,r) > 0.5);           % local indices in round r
            % [~,ord] = sort(w(idx),'descend'); idx = idx(ord); % optional
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

        %% -------- Verification (unchanged) --------
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
                kshow = min(10, size(badPairs,1));
                fprintf('  - Conflicting pairs inside rounds (showing %d/%d):\n', kshow, size(badPairs,1));
                disp(array2table(badPairs(1:kshow,:), 'VariableNames', {'Round','i','j'}));
            end
            if ~obj_ok
                fprintf('  - Objective mismatch: solver=%.6g vs recomputed=%.6g\n', objval, obj_check);
            end
        end

        fprintf('===== Run %d finished successfully. Log: %s =====\n', kRun, logPath);

    catch ME
        fprintf(2, 'ERROR in run %d: %s\n', kRun, ME.message);
        fprintf(2, 'Run parameters: startId=%d, endId=%d, p=%d, R=%d\n', startId, endId, p, R);
        fprintf(2, 'Stack:\n');
        for si = 1:numel(ME.stack)
            fprintf(2, '  %s (%d)\n', ME.stack(si).file, ME.stack(si).line);
        end
    end

    diary off;
end

%% -------- Local parser (unchanged) --------
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

