function M = run_metrics(res,varargin)
%RUN_METRICS Performance indices of one closed-loop run.
%   M = run_metrics(res,'from',10) evaluates steady-state indices for t >= 10 s.
p = inputParser;
addParameter(p,'from',[]);
addParameter(p,'tolerance',0.01);
parse(p,varargin{:});
cfg = res.cfg;
t = res.t; T = numel(t); N = size(res.x,2);
x = res.x;
if all(cfg.pin==0)
    % Leaderless: disagreement with respect to agent 1.
    e = x-x(:,1,:);
else
    e = x-reshape(res.Yref,3,1,T)-cfg.offsets;
end
position_error = squeeze(sqrt(sum(e(1:2,:,:).^2,1)));     % N x T
heading_error = squeeze(abs(e(3,:,:)));
worst = max(position_error,[],1);
from = p.Results.from;
if isempty(from)
    from = cfg.hold+(isinf(cfg.hold))*0;
    if isinf(from), from = t(end)/2; end
end
window = t >= from;
M.final_position_error = max(position_error(:,end));
M.final_heading_error = max(heading_error(:,end));
M.max_position_error_after = max(worst(window));
M.rms_position_error_after = sqrt(mean(worst(window).^2));
M.iae = cfg.Ts*sum(mean(position_error,1));
last_out = find(worst > p.Results.tolerance,1,'last');
if isempty(last_out)
    M.settling_time = 0;
elseif last_out == T
    M.settling_time = inf;
else
    M.settling_time = t(last_out+1);
end
M.control_effort = cfg.Ts*sum(res.u(:).^2)/N;
% Inter-agent distances: constrained (communicating) pairs and all pairs.
[ii,jj] = find(cfg.Ad~=0);
M.min_neighbor_distance = inf;
for q = 1:numel(ii)
    d = squeeze(sqrt(sum((x(1:2,ii(q),:)-x(1:2,jj(q),:)).^2,1)));
    M.min_neighbor_distance = min(M.min_neighbor_distance,min(d));
end
M.min_pair_distance = inf;
for a = 1:N
    for b = a+1:N
        d = squeeze(sqrt(sum((x(1:2,a,:)-x(1:2,b,:)).^2,1)));
        M.min_pair_distance = min(M.min_pair_distance,min(d));
    end
end
% Obstacle clearance margin (negative = penetration of the clearance disc).
M.min_obstacle_margin = inf;
for o = 1:numel(cfg.obstacles)
    ob = cfg.obstacles(o);
    for k = max(1,ob.appear):T
        if size(ob.center,2)==1, c = ob.center; else, c = ob.center(:,k); end
        d = sqrt(sum((x(1:2,:,k)-c).^2,1))-ob.radius;
        M.min_obstacle_margin = min(M.min_obstacle_margin,min(d));
    end
end
st = res.solve_time(:);
st = st(st>0);
if isempty(st), st = 0; end
M.solve_mean_ms = 1e3*mean(st);
M.solve_median_ms = 1e3*median(st);
M.solve_p95_ms = 1e3*prctile(st,95);
M.solve_max_ms = 1e3*max(st);
M.fallbacks = nnz(res.fallback);                 % agent-steps without a feasible hard solution
M.fallback_level2_steps = nnz(res.fallback==2);  % clearance/obstacles also relaxed
M.shifted_plan_steps = nnz(res.fallback==3);     % all fallback solutions rejected
M.nonpositive_exitflags = nnz(res.exitflag<=0 & res.solve_time>0);
if isfield(res,'violation_m')
    % Quality of the applied plans: predicted safety violation (m), residual of the
    % relaxed terminal equation and excess over the wheel-speed bound.
    M.max_predicted_violation = max(res.violation_m(:));
    M.max_predicted_compat_violation = max(res.compat_violation_m(:));
    M.max_terminal_residual = max(res.terminal_residual(:));
    M.max_input_excess = max(res.input_excess(:));
end
M.final_slack = max(vecnorm(res.slack(:,:,end),2,1));
M.position_error = position_error;
M.heading_error = heading_error;
end
