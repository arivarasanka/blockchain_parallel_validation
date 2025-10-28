%% atomated_solver_p1.m
% Run Problem 1 MILP over many [startId..endId] ranges with their own p.
% - Reads ranges from ranges.csv (columns: startId,endId,p[,slackR])
% - Writes per-iteration logs in ./logs and a summary CSV at the end.

clear; clc;

%% --- Config ---
csvFile   = 'tx_access_simple_21631019_21635079.csv';   % big input CSV
rangesCSV = 'ranges.csv';                               % list of windows to run
defaultSlackR = 5;                                      % fallback if slackR absent
summaryOut = 'summary_p1_ranges.csv';                   % summary file
logsDir    = 'logs';                                    % folder for logs
% ------------------------------

if ~exist(logsDir, 'dir'), mkdir(logsDir); end

% Read the ranges file
R = readtable(rangesCSV);
need = {'startId','endId','p'};
for k = 1:numel(need)
    if ~ismember(need{k}, R.Properties.VariableNames)
        error('ranges.csv must have column "%s"', need{k});
    end
end
hasSlack = ismember('slackR', R.Properties.VariableNames);

% Pre-read just the columns we need from the big CSV once
opts = detectImportOptions(csvFile, 'TextType','string');
needCols = ["tx_hash","blockNumber","access_read","access_write"];
opts.SelectedVariableNames = intersect(opts.VariableNames, needCols);
for c = ["tx_hash","access_read","access_write"]
    if any(strcmp(opts.VariableNames, c))
        idx = find(strcmp(opts.VariableNames, c),1);
        opts.VariableTypes{idx} = 'string';
    end
end
Tfull = readtable(csvFile, opts);
Nfile = height(Tfull);

% Summary accumulators
S_rows = [];   % will collect rows for the summary table

fprintf('Loaded main CSV with %d rows. Ranges to run: %d\n', Nfile, height(R));

for r = 1:height(R)
    startId = max(1, double(R.startId(r)));
    endId   = min(Nfile, double(R.endId(r)));
    p       = double(R.p(r));
    slackR  = defaultSlackR;
    if hasSlack && ~isnan(R.slackR(r)), slackR = double(R.slackR(r)); end

    tag = sprintf('%d_%d_p%d', startId, endId, p);
    logFile = fullfile(logsDir, sprintf('p1_%s.log', tag));
    fprintf('=== Running window [%d..%d], p=%d, slackR=%d ? log: %s ===\n', ...
            startId, endId, p, slackR, logFile);

    % Start logging for this iteration
    diary(logFile); diary on;
    try
        tStart = tic;
        % Run the window (returns metrics in struct M)
        M = run_p1_window(Tfull, startId, endId, p, slackR);
        totalTime = toc(tStart);

        fprintf('Finished window [%d..%d]: obj=%.0f, usedRounds=%d, exitflag=%d, totalTime=%.3fs\n', ...
                startId, endId, M.objval, M.usedRounds, M.exitflag, totalTime);

        % Append a summary row
        S_rows = [S_rows;
            {startId, endId, p, slackR, M.n, M.mEdges, M.L1, M.L2, M.R, M.objval, M.usedRounds, M.exitflag, totalTime, M.msg}]; 

    catch ME
        % Record failure with NaNs
        totalTime = NaN;
        fprintf(2, 'ERROR on window [%d..%d]: %s\n', startId, endId, ME.message);
        S_rows = [S_rows;
            {startId, endId, p, slackR, NaN, NaN, NaN, NaN, NaN, NaN, NaN, -999, totalTime, ME.message}]; 
    end
    diary off;
end

% Write summary CSV
Summary = cell2table(S_rows, 'VariableNames', ...
    {'startId','endId','p','slackR','n','mEdges','L1','L2','R','obj','usedRounds','exitflag','totalTime','note'});
writetable(Summary, summaryOut);
fprintf('Wrote summary to %s\n', summaryOut);

%% ---------------- Local function ----------------
function M = run_p1_window(Tfull, startId, endId, p, slackR)
    % Slice window
    T = Tfull(startId:endId, :);
    n = height(T);
    fprintf('Window [%d..%d] -> %d transactions\n', startId, endId, n);
    if n == 0, error('Empty window.'); end

    % Normalize to string columns for parsing
    if ~ismember("access_read",  T.Properties.VariableNames),  error('Column "access_read" not found.');  end
    if ~ismember("access_write", T.Properties.VariableNames),  error('Column "access_write" not found.'); end
    access_read  = string(T.access_read);
    access_write = string(T.access_write);

    parse_list = @(s) local_parse_list(s);

    % Build conflict graph (undirected): W?W or W?R (R?R allowed)
    writers = containers.Map('KeyType','char','ValueType','any');
    readers = containers.Map('KeyType','char','ValueType','any');

    for i = 1:n
        Rset = parse_list(access_read(i));
        Wset = parse_list(access_write(i));

        for t = 1:numel(Wset)
            k = char(Wset(t));
            if isKey(writers,k), writers(k) = [writers(k), i];
            else,                writers(k) = i;
            end
        end
        for t = 1:numel(Rset)
            k = char(Rset(t));
            if isKey(readers,k), readers(k) = [readers(k), i];
            else,                readers(k) = i;
            end
        end
    end

    Ei = []; Ej = [];
    wkeys = writers.keys;
    for kk = 1:numel(wkeys)
        k = wkeys{kk};
        wlist = writers(k);

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
    E_und = [Ei, Ej];
    if ~isempty(E_und), E_und = unique(E_und, 'rows'); end
    fprintf('Conflicts (undirected) built: |E| = %d\n', size(E_und,1));

    % Direct conflicts by local order i<j (acyclic precedence)
    if isempty(E_und)
        E = zeros(0,2);
    else
        E = E_und;  % already i<j
    end

    % Compute R lower bounds
    L1 = ceil(n / p);
    L2 = longest_chain_length(n, E);
    R  = max(L1, L2) + slackR;
    fprintf('Lower bounds: ceil(n/p)=%d, longest_chain=%d -> R=%d\n', L1, L2, R);

    % Solve MILP (uses your existing function)
    [xMat, y, objval, exitflag, output] = solve_min_rounds_p1(n, R, p, E); 

    usedRounds = find(sum(xMat,1) > 1e-12, 1, 'last');
    if isempty(usedRounds), usedRounds = 0; end

    % Collect metrics
    M = struct();
    M.n = n;
    M.mEdges = size(E_und,1);
    M.L1 = L1;
    M.L2 = L2;
    M.R  = R;
    M.objval = objval;
    M.usedRounds = usedRounds;
    M.exitflag = exitflag;
    if isfield(output,'message')
        M.msg = string(output.message);
    else
        M.msg = "";
    end
end

%% -------- Local: parse list of keys --------
function toks = local_parse_list(s)
    if ismissing(s) || strlength(s)==0
        toks = strings(0,1); return;
    end
    t = lower(strtrim(s));
    t = replace(t, ["[","]","{","}","'",char(34)], "");
    t = regexprep(t, '[;,\|\s]+', ' ');
    if strlength(t)==0
        toks = strings(0,1); return;
    end
    toks = string(strsplit(t, ' '));
    toks = toks(toks~="");
end
