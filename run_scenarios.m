function results = run_scenarios(names)
%RUN_SCENARIOS Simulate the three nominal scenarios of Section 5 of the paper.
%   results = run_scenarios() runs
%     'leaderless'        leaderless consensus on a directed cycle      (Section 5.2)
%     'leader_consensus'  leader-follower consensus with two obstacles  (Section 5.3)
%     'formation'         square formation following a moving leader    (Section 5.4)
%   results = run_scenarios({'formation'}) runs only the listed scenarios.
%   The closed-loop trajectories are saved in results/nominal.mat and a short
%   summary is printed; plot_results('nominal') draws the figures.

setup_paths();
if nargin < 1, names = {'leaderless','leader_consensus','formation'}; end
outdir = fullfile(fileparts(mfilename('fullpath')),'results');
if ~isfolder(outdir), mkdir(outdir); end
file = fullfile(outdir,'nominal.mat');
results = struct();
if isfile(file), S = load(file); results = S.results; end
for q = 1:numel(names)
    fprintf('Simulating %s ...\n',names{q});
    clock = tic;
    results.(names{q}) = dnmpc_sim(scenario(names{q}));
    fprintf('  done in %.1f min\n',toc(clock)/60);
end
save(file,'results','-v7.3');
fprintf('Results saved in %s\n\n',file);

fprintf('%-18s %10s %14s %12s %16s\n','Scenario','Settl. (s)','Final err (m)','Min dist (m)','Solve mean (ms)');
for q = 1:numel(names)
    M = run_metrics(results.(names{q}));
    fprintf('%-18s %10.2f %14.2e %12.3f %16.2f\n',names{q},M.settling_time, ...
        M.final_position_error,M.min_pair_distance,M.solve_mean_ms);
end
end
