function [dec,flag,iterations,elapsed] = casadi_local(P,dec0,LB,UB,level,cfg)
%CASADI_LOCAL Local OCP of one agent solved with CasADi/IPOPT (exact Hessian).
%   Same decision variables, prediction model, constraints, bounds and warm start
%   as the fmincon implementation in dnmpc_sim. IPOPT needs twice differentiable
%   functions, so the unsquared weighted norms of the stage cost are replaced by
%   the smooth surrogate sqrt(e'Me+eps^2)-eps with eps = cfg.casadi_eps.
%   flag: 1 converged, 0 iteration limit, -2 infeasibility detected, -1 other.
persistent cache
if isempty(cache), cache = containers.Map(); end
nq = size(P.received,3); nobs = numel(P.obstacles); ns = size(P.sensed,2)*P.collision;
W = P.weights;
key = sprintf('%d_%d_%d_%d_%d_%d_%d_%d_%g_%g_%g_%g_%g_%g_%g_%g_%g_%g_%s_%s', ...
    P.Np,nq,nobs,ns,P.pinned,P.collision,P.compatibility,level, ...
    W.Q(1),W.R(1),W.F(1),W.G(1),P.slack_W(1),P.Rs,P.delta_c,cfg.casadi_eps, ...
    P.penalty,P.penalty_compat,P.input_cost,P.input_norm);
if ~isKey(cache,key)
    cache(key) = build_solver(P,level,cfg,nq,nobs,ns);
end
S = cache(key);
p = pack_parameters(P,nq,nobs,ns);
clock = tic;
sol = S.solver('x0',dec0,'lbx',LB,'ubx',UB,'lbg',S.lbg,'ubg',S.ubg,'p',p);
elapsed = toc(clock);
stats = S.solver.stats();
% IPOPT relaxes the bounds by a relative 1e-8; project back onto the box.
dec = min(max(full(sol.x),LB),UB);
iterations = double(stats.iter_count);
if stats.success
    flag = 1;
elseif contains(stats.return_status,'Infeasible')
    flag = -2;
elseif contains(stats.return_status,'Maximum_Iterations')
    flag = 0;
else
    flag = -1;
end
end

function S = build_solver(P,level,cfg,nq,nobs,ns)
import casadi.*
Np = P.Np; nd = 3*Np+3;
x = SX.sym('x',nd);
p = SX.sym('p',parameter_count(Np,nq,nobs,ns));
[xi,self,received,offsets,U0,reference,leader_offset,target,centers,radii,sensed] = ...
    unpack(p,Np,nq,nobs,ns);
U = reshape(x(1:3*Np),3,Np);
xs = x(3*Np+(1:3));
% Prediction with the discrete-time model of the semi-wheel agent.
Z = SX.zeros(3,Np+1); Z(:,1) = xi;
for k = 1:Np
    th = Z(3,k);
    J = [2*cos(th), -(2*sqrt(3)/3)*cos(th+pi/6), -(2*sqrt(3)/3)*sin(th+pi/3); ...
         2*sin(th), -(2*sqrt(3)/3)*cos(th-pi/3),  (2*sqrt(3)/3)*cos(th+pi/3); ...
         1/P.R, -1/P.R, -1/P.R];
    Z(:,k+1) = Z(:,k)+P.Ts*P.r*J*U(:,k);
end
W = P.weights; e = cfg.casadi_eps;
nrm = @(v,M) sqrt(v.'*M*v+e^2)-e;
cost = xs.'*P.slack_W*xs;
for k = 1:Np                                    % column k is prediction step k-1
    cost = cost+nrm(Z(:,k)-self(:,k),W.F);
    if strcmp(P.input_cost,'deviation'), du = U(:,k)-U0(:,k); else, du = U(:,k); end
    if strcmp(P.input_norm,'norm'), cost = cost+nrm(du,W.R); else, cost = cost+du.'*W.R*du; end
    for q = 1:nq
        cost = cost+nrm(Z(:,k)-received{q}(:,k)-offsets(:,q),W.G);
    end
    if P.pinned
        cost = cost+nrm(Z(:,k)-reference(:,k)-leader_offset,W.Q);
    end
end
% Normalized safety inequalities, as in dnmpc_sim/safety_inequalities;
% compatibility rows are collected separately for the fallback weights.
c = SX.zeros(0,1); ccomp = SX.zeros(0,1);
for o = 1:nobs
    for k = 2:Np+1
        d = Z(1:2,k)-centers{o}(:,k);
        c = [c; 1-(d.'*d)/radii(o)^2]; %#ok<AGROW>
    end
end
if P.collision
    for q = 1:nq
        for k = 2:Np
            h = 2*P.Rs+(k-1)*P.delta_c; d = Z(1:2,k)-received{q}(1:2,k);
            c = [c; 1-(d.'*d)/h^2]; %#ok<AGROW>
        end
    end
    for q = 1:ns
        for k = 2:Np
            h = 2*P.Rs+P.vmax*(k-1)*P.Ts; d = Z(1:2,k)-sensed(:,q);
            c = [c; 1-(d.'*d)/h^2]; %#ok<AGROW>
        end
    end
end
if P.compatibility
    for k = 2:Np
        d = Z(1:2,k)-self(1:2,k);
        ccomp = [ccomp; (d.'*d)/P.delta_c^2-1]; %#ok<AGROW>
    end
end
g = Z(:,Np+1)-target-xs;
switch level                         % fallback levels as in dnmpc_sim/local_cost
    case 0
        hard = [c; ccomp];
    case 1
        cost = cost+P.penalty*sum1(fmax(ccomp,0).^2);
        hard = c;
    otherwise
        cost = cost+P.penalty*sum1(fmax(c,0).^2)+P.penalty_compat*sum1(fmax(ccomp,0).^2);
        hard = SX.zeros(0,1);
end
g = [g; hard];
S.lbg = [zeros(3,1); -inf(numel(hard),1)]; S.ubg = zeros(3+numel(hard),1);
nlp = struct('x',x,'p',p,'f',cost,'g',g);
ipopt = struct('print_level',0,'max_iter',cfg.max_iterations,'tol',1e-6, ...
    'constr_viol_tol',0.1*cfg.constraint_tolerance,'acceptable_tol',1e-4, ...
    'mu_init',1e-3,'sb','yes');
S.solver = nlpsol('local','ipopt',nlp,struct('print_time',false,'ipopt',ipopt));
end

function n = parameter_count(Np,nq,nobs,ns)
n = 3+3*(Np+1)+3*(Np+1)*nq+3*nq+3*Np+3*(Np+1)+3+3+2*(Np+1)*nobs+nobs+2*ns;
end

function [xi,self,received,offsets,U0,reference,leader_offset,target,centers,radii,sensed] = unpack(p,Np,nq,nobs,ns)
at = 0;
    function v = take(r,c)
        v = reshape(p(at+(1:r*c)),r,c); at = at+r*c;
    end
xi = take(3,1);
self = take(3,Np+1);
received = cell(1,nq);
for q = 1:nq, received{q} = take(3,Np+1); end
offsets = take(3,nq);
U0 = take(3,Np);
reference = take(3,Np+1);
leader_offset = take(3,1);
target = take(3,1);
centers = cell(1,nobs);
for o = 1:nobs, centers{o} = take(2,Np+1); end
radii = take(nobs,1);
sensed = take(2,ns);
end

function p = pack_parameters(P,nq,nobs,ns)
centers = zeros(0,1); radii = zeros(nobs,1);
for o = 1:nobs
    centers = [centers; P.obstacles(o).centers(:)]; %#ok<AGROW>
    radii(o) = P.obstacles(o).radius;
end
sensed = P.sensed(:,1:ns);
p = [P.xi(:); P.self(:); reshape(P.received,[],1); reshape(P.offsets(:,1:nq),[],1); P.U0(:); ...
     P.reference(:); P.leader_offset(:); P.target(:); centers; radii; sensed(:)];
end
