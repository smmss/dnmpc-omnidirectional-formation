function cfg = scenario(name,varargin)
%SCENARIO Configuration of the simulation studies reported in Section 5.
%   cfg = scenario('formation')         displacement-based formation (Sec. 5.4)
%   cfg = scenario('leader_consensus')  leader-follower consensus with obstacles (Sec. 5.3)
%   cfg = scenario('leaderless')        leaderless consensus on a directed cycle (Sec. 5.2)
%   cfg = scenario('ring','N',8)        ring formation of N agents (scalability study)
%   Name-value pairs override any field of the returned structure.

cfg = struct();
cfg.Ts = 0.01; cfg.Np = 10; cfg.R = 0.2; cfg.r = 0.04; cfg.u_max = 20;
cfg.weights = struct('Q',10*eye(3),'R',0.1*eye(3),'F',20*eye(3),'G',4*eye(3));
cfg.input_norm = 'quadratic';                               % ||u-u_ff||_R^2
speed = 0.4; radius = 2;
switch name
    case 'formation'
        % Side of the desired square. The gap between the diagonal agents 1 and 3,
        % sqrt(2)*d, exceeds twice the tightened clearance 2*Rs+(Np-1)*delta_c = 0.58 m
        % for d > 0.82 m, so that an agent can pass between them.
        p = inputParser; p.KeepUnmatched = true; addParameter(p,'d',1.2); parse(p,varargin{:});
        d = p.Results.d;
        % 1->2, 3->2, 1->3, 1->4, 2->4, 3->4: every pair exchanges trajectories
        % and the graph is acyclic.
        cfg.Ad = [0 0 0 0; 1 0 1 0; 1 0 0 0; 1 1 1 0];
        cfg.weights.F = 40*eye(3);                         % |O_1| = 3: F_1 > 36 I_3 by (67)
        cfg.pin = [1 0 1 0];
        addParameter(p,'graph','full'); parse(p,varargin{:});
        if strcmp(p.Results.graph,'chain')
            % Alternative graph for comparison: leader -> 1 -> 2 -> 3 -> 4.
            cfg.Ad = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0];
            cfg.weights.F = 20*eye(3);                     % |O_i| <= 1
            cfg.pin = [1 0 0 0];
        elseif strcmp(p.Results.graph,'chain14')
            % Chain with the redundant link 1 -> 4; |O_1| = 2: F_1 > 16 I_3.
            cfg.Ad = [0 0 0 0; 1 0 0 0; 0 1 0 0; 1 0 1 0];
            cfg.weights.F = 20*eye(3);
            cfg.pin = [1 0 0 0];
        end
        cfg.offsets = [-0.5*d 0.5*d 0.5*d -0.5*d; 0.5*d 0.5*d -0.5*d -0.5*d; 0 0 0 0];
        cfg.x0 = [-0.14 1.26 1.83 -1.05; -1.02 -0.97 -1.91 -1.95; pi/8 pi/4 pi/3 pi/6];
        cfg.hold = 2; cfg.Tf = 34; cfg.theta_ref = pi/4; cfg.center_phase = 'cos';
        cfg.slack_W = 1e5*eye(3);
        cfg.collision = struct('enabled',true,'Rs',0.2,'delta_c',0.02,'compatibility',true);
    case 'leader_consensus'
        cfg.Ad = [0 0 0 0; 1 0 0 0; 0 1 0 0; 0 0 1 0];      % cascade 1->2->3->4
        cfg.pin = [1 0 0 1];
        cfg.offsets = zeros(3,4);
        cfg.x0 = [0 0 2 2; 0 2 2 0; pi/8 pi/4 pi/3 pi/6];
        cfg.hold = 0; cfg.Tf = 32; cfg.theta_ref = 0; cfg.center_phase = 'sin';
        cfg.slack_W = 1e4*eye(3);
        cfg.weights.R = 0.05*eye(3);       % faster recovery after the obstacle detours
        clearance = 0.2;
        cfg.obstacles = struct('center',{[-2.1;0],[0;-1.8]}, ...
            'radius',{0.3+clearance,0.4+clearance},'appear',{1,1});
        cfg.collision = struct('enabled',false,'Rs',0,'delta_c',0,'compatibility',false);
    case 'leaderless'
        cfg.Ad = [0 0 0 1; 1 0 0 0; 0 1 0 0; 0 0 1 0];      % cycle 1->2->3->4->1
        cfg.pin = [0 0 0 0];
        cfg.offsets = zeros(3,4);
        cfg.x0 = [0 0 2 2; 0 2 2 0; pi/8 pi/4 pi/3 pi/6];
        cfg.hold = inf; cfg.Tf = 10; cfg.theta_ref = 0; cfg.center_phase = 'none';
        cfg.slack_W = 1e4*eye(3);
        cfg.collision = struct('enabled',false,'Rs',0,'delta_c',0,'compatibility',false);
    case 'ring'
        p = inputParser; p.KeepUnmatched = true; addParameter(p,'N',8); parse(p,varargin{:});
        N = p.Results.N;
        spacing = 0.7;
        rho = spacing/(2*sin(pi/N));
        phi = 2*pi*(0:N-1)/N;
        cfg.offsets = [rho*cos(phi); rho*sin(phi); zeros(1,N)];
        cfg.Ad = zeros(N);
        for a = 2:N
            cfg.Ad(a,a-1) = 1;
            if a >= 3, cfg.Ad(a,a-2) = 1; end
        end
        cfg.pin = zeros(1,N); cfg.pin([1 floor(N/2)+1]) = 1;
        cfg.hold = 2; cfg.Tf = 12; cfg.theta_ref = pi/4; cfg.center_phase = 'cos';
        stream = RandStream('mt19937ar','Seed',100+N);
        start = [radius;0;0]+cfg.offsets;
        cfg.x0 = start+[0.1*(2*rand(stream,2,N)-1); 0.3*rand(stream,1,N)];
        cfg.slack_W = 1e5*eye(3);
        cfg.collision = struct('enabled',true,'Rs',0.2,'delta_c',0.02,'compatibility',true);
    otherwise
        error('scenario:Name','Unknown scenario %s.',name);
end

% Overrides.
for q = 1:2:numel(varargin)
    if ~any(strcmp(varargin{q},{'N','d','graph'}))
        cfg.(varargin{q}) = varargin{q+1};
    end
end
cfg.name = name;
cfg.Yref = reference(cfg,speed,radius);
end

function Yref = reference(cfg,speed,radius)
T = round(cfg.Tf/cfg.Ts)+1;
columns = T+cfg.Np+80;
time = (0:columns-1)*cfg.Ts;
omega = speed/radius;
switch cfg.center_phase
    case 'cos'
        phase = omega*max(time-cfg.hold,0);
        Yref = [radius*cos(phase); radius*sin(phase); cfg.theta_ref*ones(1,columns)];
    case 'sin'
        phase = omega*max(time-cfg.hold,0);
        Yref = [radius*sin(phase); radius*cos(phase); cfg.theta_ref*ones(1,columns)];
    case 'none'
        Yref = zeros(3,columns);
end
end
