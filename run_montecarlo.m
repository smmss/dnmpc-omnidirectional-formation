function rows = run_montecarlo(varargin)
%RUN_MONTECARLO Monte Carlo study of the formation scenario (Section 5.8, Table 6).
%   rows = run_montecarlo() runs 1000 random initial conditions under the seven
%   conditions of Table 6 and prints the table.
%
%   Name-value options:
%     'runs'        number of random initial conditions (default 1000)
%     'workers'     parallel workers; 0 runs serially (default 8)
%     'conditions'  cell array with a subset of the conditions (default: all)
%     'graph'       formation graph: 'full' (paper, Fig. 8, default), 'chain'
%                   (0->1->2->3->4) or 'chain14' (chain plus 1->4), as compared
%                   in the response to the reviewers
%
%   Every finished run is stored in results/montecarlo_<graph>_runs/, so an
%   interrupted study resumes where it stopped. Each run takes about 1-3 min on
%   one core; 7 x 1000 runs took about 13.5 h with 8 workers. One IPOPT worker
%   needs up to about 3.5 GB of memory.
%
%   Example (quick test): run_montecarlo('runs',10,'workers',0,'conditions',{'nominal'})

p = inputParser;
addParameter(p,'runs',1000);
addParameter(p,'workers',8);
addParameter(p,'conditions',{});
addParameter(p,'graph','full');
parse(p,varargin{:});
runs = p.Results.runs; graph = p.Results.graph;
setup_paths();

variants = conditions(graph);
names = fieldnames(variants);
if ~isempty(p.Results.conditions), names = intersect(names,p.Results.conditions,'stable'); end
root = fileparts(mfilename('fullpath'));
jobdir = fullfile(root,'results',sprintf('montecarlo_%s_runs',graph));
if ~isfolder(jobdir), mkdir(jobdir); end
ics = random_initial_conditions(runs);
jobs = numel(names)*runs;
file = @(job) fullfile(jobdir,sprintf('%s_%04d.mat',names{ceil(job/runs)},job-(ceil(job/runs)-1)*runs));
todo = find(arrayfun(@(job) ~finished(file(job)),1:jobs));
fprintf('%d conditions x %d runs, %d of %d runs left.\n',numel(names),runs,numel(todo),jobs);

workers = p.Results.workers;
if workers > 0 && isempty(gcp('nocreate')), parpool('local',workers); end
parfor (q = 1:numel(todo), workers)
    job = todo(q);
    c = ceil(job/runs); k = job-(c-1)*runs;
    over = variants.(names{c});
    cfg = scenario('formation','Tf',20,'x0',ics(:,:,k),'graph',graph);
    f = fieldnames(over);
    for z = 1:numel(f), cfg.(f{z}) = over.(f{z}); end
    cfg.seed = 1000+k;
    clock = tic;
    try
        res = dnmpc_sim(cfg);
        M = run_metrics(res,'from',10);
        row = scalar_fields(M);
        row.condition = string(names{c}); row.run = k; row.ok = true;
        row.pair_min_distance = pair_distances(res);
    catch err
        row = struct('condition',string(names{c}),'run',k,'ok',false,'message',string(err.message));
    end
    save_row(file(job),row);
    fprintf('%s run %d finished (%.0f s)\n',names{c},k,toc(clock));
end

rows = cell(jobs,1);
for job = 1:jobs
    S = load(file(job),'row'); rows{job} = S.row;
end
save(fullfile(root,'results',sprintf('montecarlo_%s.mat',graph)),'rows','ics','variants','runs');
print_table(rows,names);
end

%% ------------------------------------------------------------------------
function variants = conditions(graph)
% The seven conditions of Table 6.
variants.nominal = struct();
variants.noise = struct('noise_std',[0.002 0.002 0.005]);
variants.slip_mismatch = struct('slip_max',0.10,'r_true',0.04*0.97,'R_true',0.2*1.03);
variants.delay_loss = struct('delay',2,'loss_prob',0.2);
nosensing = struct('enabled',true,'Rs',0.2,'delta_c',0.02,'compatibility',true,'sensing',false);
switch graph
    case 'full'      % links 3->2 and 1->4 fail at t = 8 s
        variants.link_failure = struct('link_failures',[801 2 3; 801 4 1]);
        variants.link_failure_nosensing = struct('link_failures',[801 2 3; 801 4 1],'collision',nosensing);
    case 'chain'     % every link is a bridge: no failure case
    case 'chain14'   % link 3->4 fails at t = 8 s or in the transient (t = 1 s)
        variants.link_failure = struct('link_failures',[801 4 3]);
        variants.link_failure_nosensing = struct('link_failures',[801 4 3],'collision',nosensing);
        variants.link_failure_early = struct('link_failures',[101 4 3]);
        variants.link_failure_early_nosensing = struct('link_failures',[101 4 3],'collision',nosensing);
    otherwise
        error('run_montecarlo:Graph','Unknown graph %s.',graph);
end
variants.combined = struct('noise_std',[0.002 0.002 0.005],'slip_max',0.10, ...
    'r_true',0.04*0.97,'R_true',0.2*1.03,'delay',2,'loss_prob',0.2);
end

function ics = random_initial_conditions(runs)
% Positions uniform in [-1.5,2] x [-2.5,-0.5] m with pairwise distances of at
% least 0.6 m, headings uniform in [0,pi/2] rad; the same for every condition.
stream = RandStream('mt19937ar','Seed',2026);
ics = zeros(3,4,runs);
for k = 1:runs
    while true
        p = [-1.5+3.5*rand(stream,1,4); -2.5+2*rand(stream,1,4)];
        D = sqrt((p(1,:)-p(1,:).').^2+(p(2,:)-p(2,:).').^2)+eye(4);
        if min(D(:)) >= 0.6, break; end
    end
    ics(:,:,k) = [p; (pi/2)*rand(stream,1,4)];
end
end

function row = scalar_fields(M)
row = struct();
f = fieldnames(M);
for q = 1:numel(f)
    if isscalar(M.(f{q})), row.(f{q}) = M.(f{q}); end
end
end

function d = pair_distances(res)
% Rows [a b smallest distance over the run, 1 if a and b exchange trajectories].
N = size(res.x,2);
pairs = nchoosek(1:N,2);
d = zeros(size(pairs,1),4);
for q = 1:size(pairs,1)
    a = pairs(q,1); b = pairs(q,2);
    gap = vecnorm(squeeze(res.x(1:2,a,:)-res.x(1:2,b,:)),2,1);
    d(q,:) = [a b min(gap) res.cfg.Ad(a,b)~=0 || res.cfg.Ad(b,a)~=0];
end
end

function ok = finished(file)
ok = false;
if isfile(file)
    S = load(file,'row');
    ok = isfield(S,'row') && S.row.ok;
end
end

function save_row(file,row)
save(file,'row');
end

function print_table(rows,names)
% Columns of Table 6. A run counts as converged if its final error is below
% 2 cm; the settling time is the first entry into the 1 cm band.
limit = 0.4*sqrt(1-1e-3)-1e-9;            % 2*Rs with the constraint tolerance
rows = rows(cellfun(@(r) r.ok,rows));
rows = [rows{:}];
fprintf(['\n%-30s %9s %7s %10s %10s %7s %7s %6s %13s %10s\n'],'Condition','Conv.','Settl.', ...
    'RMS (mm)','Final (mm)','Effort','dE (m)','d<2Rs','Solve (ms)','Fallback');
for c = 1:numel(names)
    r = rows([rows.condition]==names{c});
    r = r([r.ok]);
    if isempty(r), continue; end
    fe = [r.final_position_error];
    st = [r.settling_time];
    % Under slip the steady error is about 1 cm and most runs never enter the band.
    if nnz(isfinite(st)) < numel(st)/2, st = NaN; else, st = st(isfinite(st)); end
    D = cat(3,r.pair_min_distance); d = squeeze(D(:,3,:));
    fprintf('%-30s %4d/%-4d %7.2f %10.3g %10.3g %7.0f %7.3f %6d %6.1f/%-6.1f %5d (%d)\n', ...
        names{c},nnz(fe<0.02),numel(r),median(st),1e3*median([r.rms_position_error_after]), ...
        1e3*median(fe),median([r.control_effort]),min(d(:)),nnz(any(d<limit,1)), ...
        median([r.solve_mean_ms]),median([r.solve_p95_ms]),nnz([r.fallbacks]>0), ...
        nnz([r.fallback_level2_steps]>0));
end
end
