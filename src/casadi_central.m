function [dec,flag,iterations,elapsed] = casadi_central(P,dec0,LB,UB,cfg)
%CASADI_CENTRAL Centralized NMPC benchmark solved with CasADi/IPOPT (exact Hessian).
%   Same decision variables, cost, terminal equation and constraints as
%   central_cost/central_constraints in dnmpc_sim, with the weighted norms
%   replaced by the smooth surrogate sqrt(e'Me+eps^2)-eps as in casadi_local.
%   flag: 1 converged, 0 iteration limit, -2 infeasibility detected, -1 other.
persistent cache
if isempty(cache), cache = containers.Map(); end
N = P.N; nobs = numel(P.obstacles); W = P.weights;
key = sprintf('%d_%d_%d_%s_%s_%d_%g_%g_%g_%g_%g_%g_%g_%s_%s', ...
    N,P.Np,nobs,mat2str(P.Ad~=0),mat2str(P.pin~=0),P.collision,P.Rs, ...
    W.Q(1),W.R(1),W.G(1),P.slack_W(1),cfg.casadi_eps,cfg.max_iterations,P.input_cost,P.input_norm);
if ~isKey(cache,key)
    cache(key) = build_solver(P,cfg,nobs);
end
S = cache(key);
centers = zeros(0,1); radii = zeros(nobs,1);
for o = 1:nobs
    centers = [centers; P.obstacles(o).centers(:)]; %#ok<AGROW>
    radii(o) = P.obstacles(o).radius;
end
p = [P.xi(:); P.reference(:); P.U0(:); P.offsets(:); centers; radii];
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

function S = build_solver(P,cfg,nobs)
import casadi.*
N = P.N; Np = P.Np; nu = 3*Np*N;
x = SX.sym('x',nu+3*N);
p = SX.sym('p',3*N+3*(Np+1)+3*Np+3*N+2*(Np+1)*nobs+nobs);
at = 0;
    function v = take(r,c)
        v = reshape(p(at+(1:r*c)),r,c); at = at+r*c;
    end
xi = take(3,N); reference = take(3,Np+1); U0 = take(3,Np); offsets = take(3,N);
centers = cell(1,nobs);
for o = 1:nobs, centers{o} = take(2,Np+1); end
radii = take(nobs,1);
W = P.weights; e = cfg.casadi_eps;
nrm = @(v,M) sqrt(v.'*M*v+e^2)-e;
Z = cell(1,N); U = cell(1,N);
for a = 1:N
    U{a} = reshape(x((a-1)*3*Np+(1:3*Np)),3,Np);
    Z{a} = SX.zeros(3,Np+1); Z{a}(:,1) = xi(:,a);
    for k = 1:Np
        th = Z{a}(3,k);
        J = [2*cos(th), -(2*sqrt(3)/3)*cos(th+pi/6), -(2*sqrt(3)/3)*sin(th+pi/3); ...
             2*sin(th), -(2*sqrt(3)/3)*cos(th-pi/3),  (2*sqrt(3)/3)*cos(th+pi/3); ...
             1/P.R, -1/P.R, -1/P.R];
        Z{a}(:,k+1) = Z{a}(:,k)+P.Ts*P.r*J*U{a}(:,k);
    end
end
cost = 0; g = SX.zeros(0,1); c = SX.zeros(0,1);
for a = 1:N
    neighbors = find(P.Ad(a,:)~=0);
    for k = 1:Np
        if P.pin(a)
            cost = cost+nrm(Z{a}(:,k)-reference(:,k)-offsets(:,a),W.Q);
        end
        for j = neighbors
            cost = cost+nrm(Z{a}(:,k)-Z{j}(:,k)-(offsets(:,a)-offsets(:,j)),W.G);
        end
        if strcmp(P.input_cost,'deviation'), du = U{a}(:,k)-U0(:,k); else, du = U{a}(:,k); end
        if strcmp(P.input_norm,'norm'), cost = cost+nrm(du,W.R); else, cost = cost+du.'*W.R*du; end
    end
    s = x(nu+3*(a-1)+(1:3));
    cost = cost+s.'*P.slack_W*s;
    g = [g; Z{a}(:,end)-reference(:,end)-offsets(:,a)-s]; %#ok<AGROW>
    for o = 1:nobs
        for k = 2:Np+1
            d = Z{a}(1:2,k)-centers{o}(:,k);
            c = [c; 1-(d.'*d)/radii(o)^2]; %#ok<AGROW>
        end
    end
    if P.collision
        for j = neighbors
            for k = 2:Np
                d = Z{a}(1:2,k)-Z{j}(1:2,k);
                c = [c; 1-(d.'*d)/(2*P.Rs)^2]; %#ok<AGROW>
            end
        end
    end
end
S.lbg = [zeros(3*N,1); -inf(numel(c),1)]; S.ubg = zeros(3*N+numel(c),1);
nlp = struct('x',x,'p',p,'f',cost,'g',[g; c]);
ipopt = struct('print_level',0,'max_iter',cfg.max_iterations,'tol',1e-6, ...
    'constr_viol_tol',0.1*cfg.constraint_tolerance,'acceptable_tol',1e-4, ...
    'mu_init',1e-3,'sb','yes');
S.solver = nlpsol('central','ipopt',nlp,struct('print_time',false,'ipopt',ipopt));
end
