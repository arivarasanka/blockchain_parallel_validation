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

%% heuristic_p2_complex.m 
clear; clc;

% ---- Config ----
p        = 1;                                % #cores
R        = 630000;                         % global horizon; [] => sum(t)
txsCsv   = 'tx_access_21631019_21633879.csv';
rowFrom  = 4037; rowTo = 4280;

%% --- Inspect columns to choose the minimal set ---
opts0 = detectImportOptions(txsCsv, 'TextType','string');
vars  = string(opts0.VariableNames);

need = ["tx_hash","gas_used","access_read","access_write"];
wCol = "weight_theoretical_eth";
if ~ismember(wCol, vars)
    if ismember("weight_wei", vars)
        wCol = "weight_wei";
    elseif ismember("weight_eth", vars)
        wCol = "weight_eth";
    else
        warning('No weight column found; will default to ones.');
        wCol = "";  % handle later
    end
end
sel = need;
if strlength(wCol)>0, sel = [sel, wCol]; end

% Set only the needed columns & types
opts = opts0;
opts.SelectedVariableNames = intersect(vars, sel);
for c = ["tx_hash","access_read","access_write"]
    if any(strcmp(opts.VariableNames, c)), opts = setvartype(opts, c, 'string'); end
end

Tfull = readtable(txsCsv, opts);

% Clamp window
Nfile = height(Tfull);
rowFrom = max(1,rowFrom); rowTo = min(Nfile,rowTo);
if rowFrom>rowTo, error('Empty range'); end
T = Tfull(rowFrom:rowTo,:);
n = height(T);
fprintf('Loaded rows [%d..%d], n=%d (cols: %s)\n', rowFrom,rowTo,n, strjoin(T.Properties.VariableNames,','));

%% --- Durations t_i from gas_used ---
assert(ismember("gas_used", string(T.Properties.VariableNames)), 'gas_used missing');
t_raw = T.gas_used;
t = double(t_raw); if isstring(t_raw)||iscell(t_raw), t = str2double(string(t_raw)); end
t(~isfinite(t)) = 0; t = max(t(:),0);

%% --- Weights w_i (minimal, with fallback only if needed) ---
if strlength(wCol)==0
    w = ones(n,1);
else
    wr = T.(wCol);
    if strcmp(wCol,"weight_eth"), w = double(wr)*1e18;
    else, w = double(wr); if isstring(wr)||iscell(wr), w = str2double(string(wr)); end
    end
    w(~isfinite(w)) = 0; w = max(w(:),0);
end

if isempty(R), R = sum(t); end
fprintf('p=%d, horizon R=%g\n', p, R);
tstart=tic;

%% --- Build conflicts (only from access_read/write) ---
assert(all(ismember(["access_read","access_write"], string(T.Properties.VariableNames))), ...
       'access_read/write required');
Rcol = string(T.access_read); Wcol = string(T.access_write);

writers = containers.Map('KeyType','char','ValueType','any');
readers = containers.Map('KeyType','char','ValueType','any');

for i = 1:n
    Rset = parse_list(Rcol(i));
    Wset = parse_list(Wcol(i));
    for u = 1:numel(Wset)
        k = char(Wset(u)); writers(k) = [getdef(writers,k,[]), i];
    end
    for u = 1:numel(Rset)
        k = char(Rset(u)); readers(k) = [getdef(readers,k,[]), i];
    end
end

Ei = []; Ej = [];
wkeys = writers.keys;
for kk = 1:numel(wkeys)
    k = wkeys{kk};
    wlist = writers(k);
    % writer-writer
    m = numel(wlist);
    for a=1:max(0,m-1)
        ia = wlist(a);
        for b=a+1:m
            ib = wlist(b);
            i1=min(ia,ib); i2=max(ia,ib);
            Ei(end+1,1)=i1; Ej(end+1,1)=i2; 
        end
    end
    % writer-reader
    if isKey(readers,k)
        rlist = readers(k);
        for a=1:m
            ia = wlist(a);
            for b=1:numel(rlist)
                ib = rlist(b); if ia==ib, continue; end
                i1=min(ia,ib); i2=max(ia,ib);
                Ei(end+1,1)=i1; Ej(end+1,1)=i2; 
            end
        end
    end
end
E = [Ei, Ej]; if ~isempty(E), E = unique(E,'rows'); end
neighbors = cell(n,1);
for e=1:size(E,1), a=E(e,1); b=E(e,2); neighbors{a}(end+1)=b; neighbors{b}(end+1)=a; end
deg = cellfun(@numel, neighbors);

%% --- Static scoring & key ---
epsT = 1e-12;
rho = w ./ max(t,epsT);
s   = rho ./ (deg + 1);
ids = (1:n).';
Key = [ -s, -rho, -w, +deg, +ids ];
[~, ord] = sortrows(Key, [1 2 3 4 5]); order = ord(:);

%% --- Event-driven scheduling (simple, no deletions) ---
status      = zeros(n,1);        % 0/1/2 = UNSCHED/RUNNING/DONE
assign_core = zeros(n,1);
start_t     = zeros(n,1);
end_t       = zeros(n,1);
running     = false(n,1);
cores_busy  = false(p,1);
tcur        = 0;

% initial fill
for c=1:p
    i = next_feasible(order,status,running,neighbors,t,R,tcur);
    if i==0, break; end
    start_tx(i,c);
end

while true
    if ~any(status==1)
        if tcur>=R, break; end
        i = next_feasible(order,status,running,neighbors,t,R,tcur);
        if i==0, break; end
        c = find(~cores_busy,1);
        start_tx(i,c);
        continue;
    end
    run = find(status==1);
    [t_next,k] = min(end_t(run));
    if t_next >= R, tcur = R; break; end
    tcur = t_next;
    finset = run(end_t(run)==t_next);
    for z=1:numel(finset)
        j = finset(z);
        status(j)=2; running(j)=false; cores_busy(assign_core(j))=false;
    end
    free = find(~cores_busy);
    for u=1:numel(free)
        c = free(u);
        i = next_feasible(order,status,running,neighbors,t,R,tcur);
        if i==0, break; end
        start_tx(i,c);
    end
end

% results
sel = find(status>0);
fprintf('Scheduled %d/%d | total weight = %.8f | makespan = %.6g (<= R=%g)\n', ...
    numel(sel), n, sum(w(sel)), max([0;end_t(sel)]), R);

totalSec = toc(tstart);
fprintf('Total time: %.3f s\n', totalSec);

% show
[~, so] = sort(start_t(sel));
ii = sel(so);
Out = table(T.tx_hash(ii), start_t(ii), end_t(ii), assign_core(ii), ...
            t(ii), w(ii), deg(ii), rho(ii), s(ii), ...
            'VariableNames', {'tx_hash','start','finish','core','t','w','deg','rho','score'});
disp(Out);
disp(sum(t(ii)));

%% --- helpers ---
function val = getdef(map,k,def)
    if isKey(map,k), val = map(k); else, val = def; end
end
function toks = parse_list(s)
    if ismissing(s) || strlength(s)==0, toks = strings(0,1); return; end
    t = lower(strtrim(s));
    t = replace(t, ["[","]","{","}","'",char(34)], "");
    t = regexprep(t, '[;,\|\s]+',' ');
    if strlength(t)==0, toks = strings(0,1); return; end
    toks = string(strsplit(t,' ')); toks = toks(toks~="");
end
function i = next_feasible(order,status,running,neighbors,t,R,tcur)
    i = 0;
    for k=1:numel(order)
        ii = order(k);
        if status(ii)~=0, continue; end
        if tcur + t(ii) > R, continue; end
        nei = neighbors{ii};
        if ~isempty(nei) && any(running(nei)), continue; end
        i = ii; return;
    end
end
function start_tx(i,c)
    status      = evalin('caller','status');        status(i)=1; assignin('caller','status',status);
    assign_core = evalin('caller','assign_core');   assign_core(i)=c; assignin('caller','assign_core',assign_core);
    tcur        = evalin('caller','tcur');
    start_t     = evalin('caller','start_t');       start_t(i)=tcur; assignin('caller','start_t',start_t);
    t           = evalin('caller','t');             ti=t(i);
    end_t       = evalin('caller','end_t');         end_t(i)=tcur+ti; assignin('caller','end_t',end_t);
    running     = evalin('caller','running');       running(i)=true; assignin('caller','running',running);
    cores_busy  = evalin('caller','cores_busy');    cores_busy(c)=true; assignin('caller','cores_busy',cores_busy);
end
%% -------- Verification: edge & core constraints --------
tol = 1e-9;

% 0) Horizon & assignment sanity
bad_horizon = find(start_t < -tol | end_t > R + tol);
bad_assign  = find(assign_core < 0 | assign_core > p);  % 0 means "not scheduled", ok

% 1) No overlap on the same core
same_core_viol = [];  % rows: [core, i, j, end_i, start_j]
for c = 1:p
    idx = find(assign_core == c);
    if numel(idx) <= 1, continue; end
    [~,o] = sort(start_t(idx));
    idx = idx(o);
    for k = 2:numel(idx)
        i = idx(k-1); j = idx(k);
        if end_t(i) > start_t(j) + tol
            same_core_viol(end+1,:) = [c, i, j, end_t(i), start_t(j)]; 
        end
    end
end

% 2) Conflict edges do not overlap in time (regardless of core)
conflict_viol = [];   % rows: [i, j, start_i, end_i, start_j, end_j, core_i, core_j]
if ~isempty(E)
    for e = 1:size(E,1)
        i = E(e,1); j = E(e,2);
        % Consider only if both were scheduled
        if assign_core(i) > 0 && assign_core(j) > 0
            overlap = (start_t(i) < end_t(j) - tol) && (start_t(j) < end_t(i) - tol);
            if overlap
                conflict_viol(end+1,:) = [i, j, start_t(i), end_t(i), start_t(j), end_t(j), assign_core(i), assign_core(j)];
            end
        end
    end
end

% 3) Report
ok = isempty(bad_horizon) && isempty(bad_assign) && isempty(same_core_viol) && isempty(conflict_viol);
if ok
    fprintf('Verification: OK | horizon, per-core sequencing, and conflict edges all satisfied.\n');
else
    fprintf('Verification: issues detected \n');
    if ~isempty(bad_horizon)
        fprintf('  - Horizon violations (start<0 or end>R): %d\n', numel(bad_horizon));
    end
    if ~isempty(bad_assign)
        fprintf('  - Invalid core assignment values: %d\n', numel(bad_assign));
    end
    if ~isempty(same_core_viol)
        k = min(10, size(same_core_viol,1));
        fprintf('  - Same-core overlaps (showing %d/%d):\n', k, size(same_core_viol,1));
        disp(array2table(same_core_viol(1:k,:), 'VariableNames', ...
            {'core','i','j','end_i','start_j'}));
    end
    if ~isempty(conflict_viol)
        k = min(10, size(conflict_viol,1));
        fprintf('  - Conflicting pairs overlapped (showing %d/%d):\n', k, size(conflict_viol,1));
        disp(array2table(conflict_viol(1:k,:), 'VariableNames', ...
            {'i','j','s_i','e_i','s_j','e_j','core_i','core_j'}));
    end
end
