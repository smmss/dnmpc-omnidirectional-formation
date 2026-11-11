function res = dnmpc_sim(cfg)
%DNMPC_SIM Closed-loop simulation of the revised DNMPC and of the comparison controllers.
%   res = dnmpc_sim(cfg) simulates N semi-wheel agents. cfg is built by the
%   scenario_* functions; see fill_defaults below for every field.
%
%   cfg.mode = 'distributed' : proposed DNMPC (Algorithm 1), one local OCP per agent,
%                              all agents solve synchronously with frozen neighbor data.
%            = 'centralized' : one joint NMPC over all agents (benchmark).
%            = 'classical'   : distributed leader-following formation law with saturation.
%
%   The plant can differ from the prediction model (slip, r/R mismatch), the measured
%   state can be noisy, and messages can be delayed, lost, or links can fail.

cfg = fill_defaults(cfg);
if strcmp(cfg.solver,'casadi') && isempty(which('casadi.SX'))
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))),'casadi'));
end
n = 3; m = 3;
Np = cfg.Np; Ts = cfg.Ts;
N = size(cfg.x0,2);
T = round(cfg.Tf/Ts)+1;
assert(size(cfg.Yref,2) >= T+Np+cfg.delay+5,'dnmpc_sim:Reference','Reference is too short.');

% Disturbance realizations are drawn before the loop so that every controller
% sees the same noise, slip and packet-loss sequence for a given seed.
stream = RandStream('mt19937ar','Seed',cfg.seed);
noise = randn(stream,n,N,T).*reshape(cfg.noise_std(:),n,1,1);
slip = cfg.slip_max*rand(stream,m,N,T);
loss = rand(stream,N,N,T) < cfg.loss_prob;

res = struct();
res.cfg = cfg;
res.t = (0:T-1)*Ts;
res.x = zeros(n,N,T);
res.x(:,:,1) = cfg.x0;
res.u = zeros(m,N,T-1);
res.slack = zeros(n,N,T-1);
res.exitflag = zeros(N,T-1);
res.solve_time = zeros(N,T-1);
res.iterations = zeros(N,T-1);
res.fallback = zeros(N,T-1);              % 0: hard OCP, 1-2: fallback levels, 3: shifted previous plan
res.violation_m = zeros(N,T-1);           % largest predicted obstacle/clearance violation of the applied plan (m)
res.compat_violation_m = zeros(N,T-1);    % largest predicted violation of the compatibility constraint (m)
res.terminal_residual = zeros(N,T-1);     % residual of the (relaxed) terminal equation
res.input_excess = zeros(N,T-1);          % largest excess of a wheel speed over the bound
res.stage_cost = zeros(N,T-1);
res.status = 'completed';

switch cfg.mode
    case 'distributed'
        res = run_distributed(cfg,res,noise,slip,loss);
    case 'centralized'
        res = run_centralized(cfg,res,noise,slip);
    case 'classical'
        res = run_classical(cfg,res,noise,slip,loss);
    otherwise
        error('dnmpc_sim:Mode','Unknown mode %s.',cfg.mode);
end
res.Yref = cfg.Yref(:,1:T);
end

%% ------------------------------------------------------------------------
function res = run_distributed(cfg,res,noise,slip,loss)
n = 3; m = 3; Np = cfg.Np; Ts = cfg.Ts;
N = size(cfg.x0,2); T = numel(res.t);
Yref = cfg.Yref;
options = solver_options(cfg,m*Np+n);
LB = [-cfg.u_max*ones(m*Np,1); -inf(n,1)];
UB = [ cfg.u_max*ones(m*Np,1);  inf(n,1)];
if ~cfg.use_slack
    LB(end-2:end) = 0; UB(end-2:end) = 0;
end

x_assumed = repmat(reshape(cfg.x0,n,1,N),1,Np+1,1);
Uopt = zeros(m,Np,N);
paths = zeros(n,Np+1,N);
history = zeros(n,Np+1,N,T);          % history(:,:,j,s) = z_j^a(.|s)
stamp = zeros(N,N);                    % time stamp of the last message j->i
inbox = zeros(n,Np+1,N,N);

for i = 1:T-1
    history(:,:,:,i) = x_assumed;
    Ad = active_graph(cfg,i);
    U0 = feedforward(Yref,i,Np,Ts,cfg.R,cfg.r);
    Delta_tail = Yref(:,i+Np+1)-Yref(:,i+Np);
    if strcmp(cfg.tail,'zero')
        Delta_tail = zeros(n,1);
    end
    for a = 1:N
        neighbors = find(Ad(a,:)~=0);
        % Receive (possibly delayed or lost) assumed trajectories.
        for j = neighbors
            sent = i-cfg.delay;
            if stamp(a,j)==0
                stamp(a,j) = 1; inbox(:,:,a,j) = history(:,:,j,1);
            end
            if sent>=1 && ~loss(a,j,i) && sent>stamp(a,j)
                stamp(a,j) = sent; inbox(:,:,a,j) = history(:,:,j,sent);
            end
        end
        received = zeros(n,Np+1,numel(neighbors));
        for q = 1:numel(neighbors)
            j = neighbors(q);
            received(:,:,q) = age_shift(inbox(:,:,a,j),i-stamp(a,j),stamp(a,j),Yref,Np,cfg.tail);
        end
        % Agents whose link to agent a has failed are kept apart with their
        % measured positions (sensing-based clearance, see safety_inequalities).
        sensed = [];
        if cfg.collision.enabled && cfg.collision.sensing
            lost = find(cfg.Ad(a,:)~=0 & Ad(a,:)==0);
            sensed = res.x(1:2,lost,i)+noise(1:2,lost,i);
        end
        P = local_problem(cfg,a,i,neighbors,received,x_assumed(:,:,a),U0,res.x(:,a,i)+noise(:,a,i),sensed);
        if P.information_count==0
            error('dnmpc_sim:EmptyInformation','Agent %d has no neighbor or pinned leader at step %d.',a,i);
        end
        dec0 = [reshape(Uopt(:,:,a),[],1); zeros(n,1)];
        warm = predict(P.xi,Uopt(:,:,a),Ts,cfg.R,cfg.r);
        if cfg.use_slack
            dec0(end-2:end) = warm(:,end)-P.target;
        end
        [dec,flag,iterations,elapsed] = solve_local(P,dec0,LB,UB,options,0,cfg);
        res.solve_time(a,i) = elapsed;
        res.iterations(a,i) = iterations;
        res.exitflag(a,i) = flag;
        % The hard problem is accepted only if the returned point satisfies every
        % (normalized) constraint within the solver tolerance. Otherwise no feasible
        % solution was found, and the fallback relaxes the constraints in two levels:
        % level 1 keeps obstacle and clearance constraints hard and penalizes the
        % compatibility violations; level 2 penalizes all violations, with a much
        % larger weight on safety than on compatibility.
        level = 0;
        while ~acceptable(dec,P,level,cfg)
            level = level+1;
            if level > 2
                dec = dec0; break;                   % shifted previous plan
            end
            [dec,flag2,~,elapsed] = solve_local(P,dec0,LB,UB,options,level,cfg);
            res.solve_time(a,i) = res.solve_time(a,i)+elapsed;
            res.exitflag(a,i) = flag2;
        end
        res.fallback(a,i) = level;                   % 0: hard, 1-2: fallback levels, 3: shifted plan
        [~,ceq] = local_constraints(dec,P,0);
        [~,~,meters,kind] = safety_inequalities(predict(P.xi,reshape(dec(1:m*Np),m,Np),Ts,cfg.R,cfg.r),[],P,numel(dec));
        res.violation_m(a,i) = max([0;meters(kind<3)]);          % obstacles and clearance
        res.compat_violation_m(a,i) = max([0;meters(kind==3)]);  % compatibility
        res.terminal_residual(a,i) = max(abs(ceq));
        res.input_excess(a,i) = max(0,max(abs(dec(1:m*Np)))-cfg.u_max);
        Uopt(:,:,a) = reshape(dec(1:m*Np),m,Np);
        res.slack(:,a,i) = dec(end-2:end);
        res.stage_cost(a,i) = local_cost(dec,P,false)-dec(end-2:end).'*cfg.slack_W*dec(end-2:end);
        paths(:,:,a) = predict(P.xi,Uopt(:,:,a),Ts,cfg.R,cfg.r);
    end
    % Apply all first inputs together, then broadcast all shifted plans.
    for a = 1:N
        res.u(:,a,i) = Uopt(:,1,a);
        res.x(:,a,i+1) = plant_step(res.x(:,a,i),Uopt(:,1,a),slip(:,a,i),cfg);
        if strcmp(cfg.tail,'zero')
            u_tail = zeros(m,1);
        else
            u_tail = increment_input(paths(3,end,a),Delta_tail,Ts,cfg.R,cfg.r);
        end
        x_assumed(:,:,a) = [paths(:,2:end,a),paths(:,end,a)+Delta_tail];
        Uopt(:,:,a) = [Uopt(:,2:end,a),u_tail];
    end
end
end

function P = local_problem(cfg,a,i,neighbors,received,own_assumed,U0,xi,sensed)
n = 3; Np = cfg.Np;
if nargin<9, sensed = zeros(2,0); end
P.sensed = reshape(sensed,2,[]); P.vmax = cfg.vmax;
P.xi = xi; P.Np = Np; P.Ts = cfg.Ts; P.R = cfg.R; P.r = cfg.r;
P.weights = cfg.weights; P.slack_W = cfg.slack_W;
P.input_cost = cfg.input_cost; P.input_norm = cfg.input_norm; P.norm_eps = cfg.norm_eps;
P.self = own_assumed;
P.received = received;
P.offsets = cfg.offsets(:,a)-cfg.offsets(:,neighbors);
P.U0 = U0;
P.pinned = cfg.pin(a)~=0;
P.reference = cfg.Yref(:,i:i+Np);
P.leader_offset = cfg.offsets(:,a);
P.information_count = numel(neighbors)+cfg.pin(a);
target = zeros(n,1);
for q = 1:numel(neighbors)
    target = target+received(:,end,q)+P.offsets(:,q);
end
if P.pinned
    target = target+P.reference(:,end)+P.leader_offset;
end
P.target = target/max(P.information_count,1);
P.obstacles = active_obstacles(cfg,i);
P.collision = cfg.collision.enabled;
P.Rs = cfg.collision.Rs; P.delta_c = cfg.collision.delta_c;
P.compatibility = cfg.collision.enabled && cfg.collision.compatibility;
P.penalty = cfg.soft_penalty; P.penalty_compat = cfg.soft_penalty_compat;
end

function [J,dJ] = local_cost(dec,P,level)
% level 0: hard problem; 1: compatibility penalized; 2: all safety inequalities penalized.
Np = P.Np;
U = reshape(dec(1:3*Np),3,Np);
xi_s = dec(3*Np+(1:3));
[Z,S] = predict(P.xi,U,P.Ts,P.R,P.r);
W = P.weights;
J = xi_s.'*P.slack_W*xi_s;
dJ = zeros(numel(dec),1);
for k = 1:Np                         % column k is prediction step k-1
    [v,gs] = weighted_norm(Z(:,k)-P.self(:,k),W.F,P.norm_eps);
    J = J+v;
    if strcmp(P.input_cost,'deviation')
        du = U(:,k)-P.U0(:,k);
    else
        du = U(:,k);
    end
    if strcmp(P.input_norm,'norm')
        [v,gu] = weighted_norm(du,W.R,P.norm_eps);
    else
        v = du.'*W.R*du; gu = 2*W.R*du;
    end
    J = J+v;
    for q = 1:size(P.received,3)
        [v,g] = weighted_norm(Z(:,k)-P.received(:,k,q)-P.offsets(:,q),W.G,P.norm_eps);
        J = J+v; gs = gs+g;
    end
    if P.pinned
        [v,g] = weighted_norm(Z(:,k)-P.reference(:,k)-P.leader_offset,W.Q,P.norm_eps);
        J = J+v; gs = gs+g;
    end
    dJ(1:3*Np) = dJ(1:3*Np)+S(:,:,k).'*gs;
    idx = 3*(k-1)+(1:3);
    dJ(idx) = dJ(idx)+gu;
end
dJ(3*Np+(1:3)) = 2*P.slack_W*xi_s;
if level > 0
    % Fallback: quadratic penalties on the violations of the relaxed inequalities.
    % Level 1 relaxes only the compatibility constraint (weight rho); level 2 also
    % relaxes obstacles and clearance (weight rho), with compatibility weighted by
    % rho_c << rho, so that safety keeps priority over plan consistency.
    [c,dC,~,kind] = safety_inequalities(Z,S,P,numel(dec));
    if level == 1
        weight = P.penalty*(kind==3);
    else
        weight = P.penalty*ones(size(c)); weight(kind==3) = P.penalty_compat;
    end
    active = max(c,0);
    J = J+sum(weight.*active.^2);
    dJ = dJ+2*(dC*(weight.*active));
end
end

function [c,ceq,dC,dCeq] = local_constraints(dec,P,level)
% Hard inequalities of the given level (see local_cost) and the terminal equation.
Np = P.Np;
U = reshape(dec(1:3*Np),3,Np);
xi_s = dec(3*Np+(1:3));
[Z,S] = predict(P.xi,U,P.Ts,P.R,P.r);
ceq = Z(:,end)-P.target-xi_s;
dCeq = [S(:,:,end).'; -eye(3)];
if level >= 2
    c = zeros(0,1); dC = zeros(numel(dec),0);
else
    [c,dC,~,kind] = safety_inequalities(Z,S,P,numel(dec));
    if level == 1
        c = c(kind<3); dC = dC(:,kind<3);
    end
end
end

function [c,dC,meters,kind] = safety_inequalities(Z,S,P,nd)
% Obstacles (steps 1..Np), tightened inter-agent clearance, sensing-based clearance
% to agents whose link has failed, and compatibility (steps 1..Np-1). Step 0 is
% fixed by the measured state and is not a decision. Every inequality c <= 0 is
% a squared distance divided by the squared radius it is compared with, so that
% the constraint tolerance corresponds to the same relative distance error for
% all constraints. `meters` returns the violations in metres (positive = violated).
Np = P.Np;
nobs = numel(P.obstacles);
nnb = size(P.received,3)*P.collision;
nsens = size(P.sensed,2)*P.collision;
ncomp = P.compatibility;
c = zeros(nobs*Np+(nnb+nsens+ncomp)*(Np-1),1);
meters = zeros(size(c));
kind = [ones(nobs*Np,1); 2*ones((nnb+nsens)*(Np-1),1); 3*ones(ncomp*(Np-1),1)];  % obstacle, clearance, compatibility
gradients = ~isempty(S);
if gradients, dC = zeros(nd,numel(c)); else, dC = []; end
row = 0;
for o = 1:nobs
    C = P.obstacles(o).centers;          % 2 x (Np+1)
    rad = P.obstacles(o).radius;
    for k = 2:Np+1
        delta = Z(1:2,k)-C(:,k);
        row = row+1;
        c(row) = 1-(delta.'*delta)/rad^2;
        meters(row) = rad-norm(delta);
        if gradients, dC(1:3*Np,row) = -2*S(1:2,:,k).'*delta/rad^2; end
    end
end
if P.collision
    for q = 1:size(P.received,3)
        for k = 2:Np
            h = 2*P.Rs+(k-1)*P.delta_c;
            delta = Z(1:2,k)-P.received(1:2,k,q);
            row = row+1;
            c(row) = 1-(delta.'*delta)/h^2;
            meters(row) = h-norm(delta);
            if gradients, dC(1:3*Np,row) = -2*S(1:2,:,k).'*delta/h^2; end
        end
    end
    % A sensed agent can be anywhere within vmax*(k-1)*Ts of its measured position.
    for q = 1:size(P.sensed,2)
        for k = 2:Np
            h = 2*P.Rs+P.vmax*(k-1)*P.Ts;
            delta = Z(1:2,k)-P.sensed(:,q);
            row = row+1;
            c(row) = 1-(delta.'*delta)/h^2;
            meters(row) = h-norm(delta);
            if gradients, dC(1:3*Np,row) = -2*S(1:2,:,k).'*delta/h^2; end
        end
    end
end
if P.compatibility
    for k = 2:Np
        delta = Z(1:2,k)-P.self(1:2,k);
        row = row+1;
        c(row) = (delta.'*delta)/P.delta_c^2-1;
        meters(row) = norm(delta)-P.delta_c;
        if gradients, dC(1:3*Np,row) = 2*S(1:2,:,k).'*delta/P.delta_c^2; end
    end
end
end

function ok = acceptable(dec,P,level,cfg)
% A solution is applied only if it is finite, respects the input bounds and the
% relaxed terminal equation, and satisfies the hard inequalities of its level
% within the constraint tolerance.
m = 3; Np = P.Np;
ok = all(isfinite(dec)) && max(abs(dec(1:m*Np))) <= cfg.u_max*(1+1e-9);
if ~ok, return; end
[c,ceq] = local_constraints(dec,P,level);
ok = max(abs(ceq)) <= cfg.constraint_tolerance && ...
    (isempty(c) || max(c) <= cfg.constraint_tolerance);
end

function [dec,flag,iterations,elapsed] = solve_local(P,dec0,LB,UB,options,level,cfg)
% Local OCP of the given level with fmincon (SQP, analytic gradients) or CasADi/IPOPT.
if strcmp(cfg.solver,'casadi')
    [dec,flag,iterations,elapsed] = casadi_local(P,dec0,LB,UB,level,cfg);
    return
end
FCN = @(d) local_cost(d,P,level);
CNS = @(d) local_constraints(d,P,level);
clock = tic;
[dec,~,flag,out] = fmincon(FCN,dec0,[],[],[],[],LB,UB,CNS,options);
elapsed = toc(clock);
iterations = out.iterations;
end

%% ------------------------------------------------------------------------
function res = run_centralized(cfg,res,noise,slip)
n = 3; m = 3; Np = cfg.Np; Ts = cfg.Ts;
N = size(cfg.x0,2); T = numel(res.t);
Yref = cfg.Yref;
nu = m*Np*N;
options = solver_options(cfg,nu+n*N);
LB = [-cfg.u_max*ones(nu,1); -inf(n*N,1)];
UB = [ cfg.u_max*ones(nu,1);  inf(n*N,1)];
if ~cfg.use_slack
    LB(nu+1:end) = 0; UB(nu+1:end) = 0;
end
U = zeros(m,Np,N);
for i = 1:T-1
    Ad = active_graph(cfg,i);
    P.Ad = Ad; P.N = N; P.Np = Np; P.Ts = Ts; P.R = cfg.R; P.r = cfg.r;
    P.weights = cfg.weights; P.slack_W = cfg.slack_W; P.pin = cfg.pin;
    P.offsets = cfg.offsets; P.reference = Yref(:,i:i+Np);
    P.U0 = feedforward(Yref,i,Np,Ts,cfg.R,cfg.r);
    P.input_cost = cfg.input_cost; P.input_norm = cfg.input_norm; P.norm_eps = cfg.norm_eps;
    P.xi = res.x(:,:,i)+noise(:,:,i);
    P.obstacles = active_obstacles(cfg,i);
    P.collision = cfg.collision.enabled; P.Rs = cfg.collision.Rs;
    dec0 = [U(:); zeros(n*N,1)];
    for a = 1:N
        warm = predict(P.xi(:,a),U(:,:,a),Ts,cfg.R,cfg.r);
        dec0(nu+3*(a-1)+(1:3)) = warm(:,end)-P.reference(:,end)-P.offsets(:,a);
    end
    if strcmp(cfg.solver,'casadi')
        [dec,flag,iterations,elapsed] = casadi_central(P,dec0,LB,UB,cfg);
    else
        FCN = @(d) central_cost(d,P);
        CNS = @(d) central_constraints(d,P);
        clock = tic;
        [dec,~,flag,out] = fmincon(FCN,dec0,[],[],[],[],LB,UB,CNS,options);
        elapsed = toc(clock); iterations = out.iterations;
    end
    res.solve_time(:,i) = elapsed;
    res.iterations(:,i) = iterations;
    res.exitflag(:,i) = flag;
    U = reshape(dec(1:nu),m,Np,N);
    res.slack(:,:,i) = reshape(dec(nu+1:end),n,N);
    for a = 1:N
        res.u(:,a,i) = U(:,1,a);
        res.x(:,a,i+1) = plant_step(res.x(:,a,i),U(:,1,a),slip(:,a,i),cfg);
        tail = increment_input(Yref(3,i+Np),Yref(:,i+Np+1)-Yref(:,i+Np),Ts,cfg.R,cfg.r);
        U(:,:,a) = [U(:,2:end,a),tail];
    end
end
end

function [J,dJ] = central_cost(dec,P)
N = P.N; Np = P.Np; nu = 3*Np*N; W = P.weights;
U = reshape(dec(1:nu),3,Np,N);
Z = zeros(3,Np+1,N); S = zeros(3,3*Np,Np+1,N);
for a = 1:N
    [Z(:,:,a),S(:,:,:,a)] = predict(P.xi(:,a),U(:,:,a),P.Ts,P.R,P.r);
end
J = 0; dJ = zeros(numel(dec),1);
for a = 1:N
    ia = (a-1)*3*Np+(1:3*Np);
    neighbors = find(P.Ad(a,:)~=0);
    for k = 1:Np
        gs = zeros(3,1);
        if P.pin(a)
            [v,g] = weighted_norm(Z(:,k,a)-P.reference(:,k)-P.offsets(:,a),W.Q,P.norm_eps);
            J = J+v; gs = gs+g;
        end
        for j = neighbors
            [v,g] = weighted_norm(Z(:,k,a)-Z(:,k,j)-(P.offsets(:,a)-P.offsets(:,j)),W.G,P.norm_eps);
            J = J+v; gs = gs+g;
            ij = (j-1)*3*Np+(1:3*Np);
            dJ(ij) = dJ(ij)-S(:,:,k,j).'*g;
        end
        dJ(ia) = dJ(ia)+S(:,:,k,a).'*gs;
        if strcmp(P.input_cost,'deviation')
            du = U(:,k,a)-P.U0(:,k);
        else
            du = U(:,k,a);
        end
        if strcmp(P.input_norm,'norm')
            [v,gu] = weighted_norm(du,W.R,P.norm_eps);
        else
            v = du.'*W.R*du; gu = 2*W.R*du;
        end
        J = J+v;
        idx = (a-1)*3*Np+3*(k-1)+(1:3);
        dJ(idx) = dJ(idx)+gu;
    end
    s = dec(nu+3*(a-1)+(1:3));
    J = J+s.'*P.slack_W*s;
    dJ(nu+3*(a-1)+(1:3)) = 2*P.slack_W*s;
end
end

function [c,ceq,dC,dCeq] = central_constraints(dec,P)
N = P.N; Np = P.Np; nu = 3*Np*N; nd = numel(dec);
U = reshape(dec(1:nu),3,Np,N);
Z = zeros(3,Np+1,N); S = zeros(3,3*Np,Np+1,N);
for a = 1:N
    [Z(:,:,a),S(:,:,:,a)] = predict(P.xi(:,a),U(:,:,a),P.Ts,P.R,P.r);
end
ceq = zeros(3*N,1); dCeq = zeros(nd,3*N);
for a = 1:N
    rows = 3*(a-1)+(1:3);
    ceq(rows) = Z(:,end,a)-P.reference(:,end)-P.offsets(:,a)-dec(nu+rows);
    dCeq((a-1)*3*Np+(1:3*Np),rows) = S(:,:,end,a).';
    dCeq(nu+rows,rows) = -eye(3);
end
c = zeros(0,1); dC = zeros(nd,0);
for a = 1:N
    ia = (a-1)*3*Np+(1:3*Np);
    for o = 1:numel(P.obstacles)
        C = P.obstacles(o).centers; rad = P.obstacles(o).radius;
        for k = 2:Np+1
            delta = Z(1:2,k,a)-C(:,k);
            c(end+1,1) = 1-(delta.'*delta)/rad^2; %#ok<AGROW>
            g = zeros(nd,1); g(ia) = -2*S(1:2,:,k,a).'*delta/rad^2;
            dC(:,end+1) = g; %#ok<AGROW>
        end
    end
    if P.collision
        for j = find(P.Ad(a,:)~=0)
            ij = (j-1)*3*Np+(1:3*Np);
            for k = 2:Np
                delta = Z(1:2,k,a)-Z(1:2,k,j);
                c(end+1,1) = 1-(delta.'*delta)/(2*P.Rs)^2; %#ok<AGROW>
                g = zeros(nd,1);
                g(ia) = -2*S(1:2,:,k,a).'*delta/(2*P.Rs)^2;
                g(ij) = 2*S(1:2,:,k,j).'*delta/(2*P.Rs)^2;
                dC(:,end+1) = g; %#ok<AGROW>
            end
        end
    end
end
end

%% ------------------------------------------------------------------------
function res = run_classical(cfg,res,noise,slip,~)
% Distributed leader-following formation law with feedforward and saturation:
% v_i = Delta_0/Ts + k ( p_i (z_0 + d_i0 - z_i) + sum_j (z_j + d_ij - z_i) ) / |I_i|.
Ts = cfg.Ts; N = size(cfg.x0,2); T = numel(res.t);
Yref = cfg.Yref;
for i = 1:T-1
    Ad = active_graph(cfg,i);
    measured = res.x(:,:,i)+noise(:,:,i);
    delayed = res.x(:,:,max(1,i-cfg.delay))+noise(:,:,max(1,i-cfg.delay));
    clock = tic;
    for a = 1:N
        neighbors = find(Ad(a,:)~=0);
        e = zeros(3,1);
        if cfg.pin(a)
            e = e+Yref(:,i)+cfg.offsets(:,a)-measured(:,a);
        end
        for j = neighbors
            e = e+delayed(:,j)+cfg.offsets(:,a)-cfg.offsets(:,j)-measured(:,a);
        end
        e = e/max(numel(neighbors)+cfg.pin(a),1);
        v = (Yref(:,i+1)-Yref(:,i))/Ts+cfg.classical_gain*e;
        u = (cfg.r*Js(measured(3,a),cfg.R))\v;
        scale = max(1,max(abs(u))/cfg.u_max);
        res.u(:,a,i) = u/scale;
        res.x(:,a,i+1) = plant_step(res.x(:,a,i),res.u(:,a,i),slip(:,a,i),cfg);
    end
    res.solve_time(:,i) = toc(clock)/N;
end
end

%% ------------------------------------------------------------------------
function cfg = fill_defaults(cfg)
defaults = struct( ...
    'mode','distributed','Ts',0.01,'Np',10,'R',0.2,'r',0.04,'u_max',20, ...
    'weights',struct('Q',10*eye(3),'R',0.01*eye(3),'F',20*eye(3),'G',4*eye(3)), ...
    'slack_W',1e5*eye(3),'use_slack',true, ...
    'input_cost','deviation','input_norm','norm','norm_eps',0,'tail','feedforward', ...
    'obstacles',struct('center',{},'radius',{},'appear',{}), ...
    'collision',struct('enabled',false,'Rs',0.2,'delta_c',0.02,'compatibility',true), ...
    'noise_std',[0 0 0],'slip_max',0,'r_true',[],'R_true',[], ...
    'delay',0,'loss_prob',0,'link_failures',zeros(0,3), ...
    'soft_penalty',1e6,'soft_penalty_compat',1,'classical_gain',2,'seed',1, ...
    'optimality_tolerance',1e-3,'constraint_tolerance',1e-3,'step_tolerance',1e-4, ...
    'max_iterations',200,'solver','casadi','casadi_eps',1e-4,'vmax',[]);
names = fieldnames(defaults);
for q = 1:numel(names)
    if ~isfield(cfg,names{q}) || isempty(cfg.(names{q}))
        cfg.(names{q}) = defaults.(names{q});
    end
end
if isempty(cfg.r_true), cfg.r_true = cfg.r; end
if isempty(cfg.R_true), cfg.R_true = cfg.R; end
if ~isfield(cfg.collision,'sensing'), cfg.collision.sensing = true; end
if isempty(cfg.vmax)
    % Largest planar speed of an agent: the input box is a polytope, so the
    % maximum of |v_xy| is attained at one of its 8 vertices (4 r u_max here).
    V = cfg.u_max*(2*(dec2bin(0:7)-'0').'-1);
    v = cfg.r*Js(0,cfg.R)*V;
    cfg.vmax = max(vecnorm(v(1:2,:),2,1));
end
end

function options = solver_options(cfg,nd)
options = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'OptimalityTolerance',cfg.optimality_tolerance,'ConstraintTolerance',cfg.constraint_tolerance, ...
    'StepTolerance',cfg.step_tolerance,'MaxIterations',cfg.max_iterations, ...
    'MaxFunctionEvaluations',200*nd, ...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
if isfield(cfg,'check_gradients') && cfg.check_gradients
    options = optimoptions(options,'CheckGradients',true,'FiniteDifferenceType','central','Display','final');
end
end

function Ad = active_graph(cfg,i)
Ad = cfg.Ad;
for q = 1:size(cfg.link_failures,1)
    if i >= cfg.link_failures(q,1)
        Ad(cfg.link_failures(q,2),cfg.link_failures(q,3)) = 0;
    end
end
end

function obstacles = active_obstacles(cfg,i)
Np = cfg.Np;
obstacles = struct('centers',{},'radius',{});
for o = 1:numel(cfg.obstacles)
    ob = cfg.obstacles(o);
    if ~isempty(ob.appear) && i < ob.appear
        continue;
    end
    if size(ob.center,2)==1
        centers = repmat(ob.center,1,Np+1);
    else
        centers = ob.center(:,i:i+Np);
    end
    obstacles(end+1) = struct('centers',centers,'radius',ob.radius); %#ok<AGROW>
end
end

function traj = age_shift(message,age,sent,Yref,Np,tail)
% Shift a message of age `age` steps to the current time; beyond its horizon the
% message is extended with the planned leader increments (or held for a zero tail).
traj = zeros(3,Np+1);
for k = 0:Np
    idx = k+age;
    if idx <= Np
        traj(:,k+1) = message(:,idx+1);
    elseif strcmp(tail,'zero')
        traj(:,k+1) = message(:,Np+1);
    else
        traj(:,k+1) = message(:,Np+1)+Yref(:,sent+idx)-Yref(:,sent+Np);
    end
end
end

function U0 = feedforward(Yref,i,Np,Ts,R,r)
U0 = zeros(3,Np);
for k = 1:Np
    U0(:,k) = increment_input(Yref(3,i+k-1),Yref(:,i+k)-Yref(:,i+k-1),Ts,R,r);
end
end

function x_next = plant_step(x,u,slip,cfg)
% Plant with optional roller slip and geometric mismatch.
x_next = x+cfg.Ts*cfg.r_true*Js(x(3),cfg.R_true)*(u.*(1-slip));
end

function [Z,S] = predict(xi,U,Ts,R,r)
H = size(U,2);
Z = zeros(3,H+1); Z(:,1) = xi;
if nargout>1
    S = zeros(3,numel(U),H+1);
end
for k = 1:H
    eta = Ts*r*Js(Z(3,k),R);
    delta = eta*U(:,k);
    Z(:,k+1) = Z(:,k)+delta;
    if nargout>1
        A = eye(3); A(1,3) = -delta(2); A(2,3) = delta(1);
        S(:,:,k+1) = A*S(:,:,k);
        idx = 3*(k-1)+(1:3);
        S(:,idx,k+1) = S(:,idx,k+1)+eta;
    end
end
end

function J = Js(theta,R)
J = [2*cos(theta) -(2*sqrt(3)/3)*cos(theta+pi/6) -(2*sqrt(3)/3)*sin(theta+pi/3);
     2*sin(theta) -(2*sqrt(3)/3)*cos(theta-pi/3)  (2*sqrt(3)/3)*cos(theta+pi/3);
     1/R -1/R -1/R];
end

function u = increment_input(theta,delta,Ts,R,r)
u = (Ts*r*Js(theta,R))\delta;
end

function [value,gradient] = weighted_norm(e,M,epsilon)
% Weighted norm sqrt(e'Me); epsilon > 0 uses the smooth surrogate
% sqrt(e'Me+epsilon^2)-epsilon, which removes the kink at e = 0.
if nargin<3, epsilon = 0; end
value = sqrt(max(0,e.'*M*e)+epsilon^2);
if value>0
    gradient = M*e/value;
else
    gradient = zeros(size(e));
end
value = value-epsilon;
end
