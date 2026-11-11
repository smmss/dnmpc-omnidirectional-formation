function plot_results(which)
%PLOT_RESULTS Figures of the nominal scenarios and of the Monte Carlo study.
%   plot_results('nominal')     figures of Sections 5.2-5.4 from results/nominal.mat
%   plot_results('montecarlo')  error and solve-time distributions from
%                               results/montecarlo_full.mat
%   plot_results('all')         both (default)
%   The figures are written as EPS files to results/figures/<scenario>/.
if nargin<1, which = 'all'; end
setup_paths();
data = fullfile(fileparts(mfilename('fullpath')),'results');
img = fullfile(data,'figures');
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');
colors = lines(7);
agents = {'Agent 1','Agent 2','Agent 3','Agent 4'};

if any(strcmp(which,{'nominal','all'}))
    S = load(fullfile(data,'nominal.mat')); R = S.results;

    % Leaderless consensus (Section 5.2).
    res = R.leaderless; d = fullfile(img,'ConsensusLeaderless');
    t = res.t; T = numel(t);
    f = newfig(16,14);
    plot(squeeze(res.x(1,:,:)).',squeeze(res.x(2,:,:)).','LineWidth',1.5);
    grid on; axis equal;
    legend(agents,'Location','best');
    xlabel('$x(t)\;(\mathrm{m})$'); ylabel('$y(t)\;(\mathrm{m})$');
    title('Agents trajectories in the $x$--$y$ plane');
    style(gca);
    export(f,d,'Fig4_xy_trajectories');
    for w = 1:3
        f = newfig(16,12);
        plot(t(1:end-1),squeeze(res.u(w,:,:)).','LineWidth',1.5);
        grid on; legend(agents,'Location','best');
        xlabel('Time (s)'); ylabel(sprintf('$u_%d(t)\\;(\\mathrm{rad/s})$',w));
        title(sprintf('Angular velocity of wheel %d for all agents vs. time',w));
        style(gca);
        export(f,d,sprintf('Fig%d_u%d_inputs',w+4,w));
    end
    % Slack and heading panels.
    series = {squeeze(vecnorm(res.slack(1:2,:,:),2,1)), squeeze(abs(res.slack(3,:,:))), ...
        squeeze(res.x(3,:,:))};
    times = {t(1:end-1),t(1:end-1),t};
    names = {'terminal_position_slack','terminal_heading_slack','heading_angles'};
    ylabels = {'$\|\xi_{i,xy}(t)\|_2\;(\mathrm{m})$','$|\xi_{i,\theta}(t)|\;(\mathrm{rad})$', ...
        '$\theta_i(t)\;(\mathrm{rad})$'};
    for k = 1:3
        f = figure('Visible','off','Color','white','Units','centimeters', ...
            'Position',[2 2 6.2 5.2],'Renderer','painters');
        ax = axes(f,'Units','normalized','Position',[0.22 0.20 0.74 0.75], ...
            'FontName','Times New Roman','FontSize',8,'LineWidth',0.6, ...
            'TickLabelInterpreter','latex','ColorOrder',colors(1:4,:),'NextPlot','add');
        plot(ax,times{k},series{k}.','LineWidth',1.0);
        grid(ax,'on'); box(ax,'on');
        xlim(ax,[0 t(end)]); xticks(ax,0:2:t(end));
        xlabel(ax,'Time (s)','FontSize',9); ylabel(ax,ylabels{k},'FontSize',9);
        legend(ax,agents,'Location','northeast','FontSize',6.5,'Box','on');
        if k<3, ylim(ax,[0 max(series{k},[],'all')*1.08]); end
        export(f,d,names{k});
    end

    % Leader-follower consensus with obstacles (Section 5.3).
    res = R.leader_consensus; d = fullfile(img,'ConsensusLeader');
    T = numel(res.t); x = res.x;
    f = newfig(16,14); hold on; grid on;
    phi = linspace(0,2*pi,100);
    ob = res.cfg.obstacles;
    clearance = 0.2;                       % the obstacle radius in cfg includes the clearance
    ring = cell(1,2);
    for o = 1:2
        rad = ob(o).radius-clearance;
        ring{o} = [ob(o).center(1)+rad*cos(phi); ob(o).center(2)+rad*sin(phi)];
        fill(ring{o}(1,:),ring{o}(2,:),[0.925 0.925 0.925],'EdgeColor','none');
    end
    h_ref = plot(res.Yref(1,1:T-1),res.Yref(2,1:T-1),'--','Color',[0.4 0.7 0.3],'LineWidth',1.5);
    h = gobjects(1,4);
    for a = 1:4
        h(a) = plot(squeeze(x(1,a,1:T-1)),squeeze(x(2,a,1:T-1)),'Color',colors(a,:),'LineWidth',1.5);
    end
    h_ob1 = plot(ring{1}(1,:),ring{1}(2,:),'--','Color',[0.2 0.6 0.8],'LineWidth',1.5);
    h_ob2 = plot(ring{2}(1,:),ring{2}(2,:),'--','Color',[0.6 0.1 0.2],'LineWidth',1.5);
    plot(squeeze(x(1,:,1)),squeeze(x(2,:,1)),'o','MarkerSize',6,'MarkerFaceColor','w', ...
        'LineWidth',1.2,'Color','k','HandleVisibility','off');
    for a = 1:4
        plot(x(1,a,T-1),x(2,a,T-1),'s','MarkerSize',7,'MarkerFaceColor',colors(a,:), ...
            'Color',colors(a,:),'HandleVisibility','off');
    end
    h_lgd = legend([h,h_ref,h_ob1,h_ob2], ...
        {'Agent 1','Agent 2','Agent 3','Agent 4','Virtual Leader','Obstacle 1','Obstacle 2'});
    set(h_lgd,'Position',[0.39 0.62 0.18 0.25]);
    xlabel('$x(t)\;(\mathrm{m})$'); ylabel('$y(t)\;(\mathrm{m})$');
    title('Agents trajectories and obstacles in the $x$-$y$ plane');
    axis equal; style(gca);
    export(f,d,'Fig3_xy_trajectories');

    % Formation (Section 5.4).
    res = R.formation; d = fullfile(img,'FormationDisp');
    t = res.t; T = numel(t); x = res.x;
    period = 2*pi/(0.4/2);
    times = [0, res.cfg.hold+period*[1/3 2/3 1]];
    names = {'Fig1_snapshot_initial','Fig2_snapshot_one_third_lap', ...
        'Fig3_snapshot_two_thirds_lap','Fig4_snapshot_full_lap'};
    for q = 1:4
        [~,k] = min(abs(t-times(q)));
        f = newfig(14,12);
        colors = get(gca,'ColorOrder');
        h = plot(res.Yref(1,1:T),res.Yref(2,1:T),'--','LineWidth',1.5,'Color',colors(1,:)); hold on;
        % Communication links j -> i (Ad(i,j) = 1) and leader links to the pinned agents.
        pos = x(1:2,:,k);
        [to,from] = find(res.cfg.Ad);
        ref = res.Yref(1:2,k);
        for e = 1:numel(to)
            % A link of a pinned agent that is collinear with its leader link
            % (e.g. 1 -> 3, since the leader is at the centre of the square)
            % is drawn as an arc, so that it does not cover the leader link.
            p = pos(:,from(e)); v = pos(:,to(e))-p;
            bend = 0;
            for a = reshape(intersect(find(res.cfg.pin),[from(e) to(e)]),1,[])
                dl = pos(:,a)-ref;
                if abs(dl(1)*v(2)-dl(2)*v(1)) < 0.2*norm(dl)*norm(v)
                    bend = 0.3*norm(v);
                end
            end
            hl = link_arrow(p,pos(:,to(e)),'-',[0.25 0.25 0.25],bend);
        end
        for a = find(res.cfg.pin)
            hp = link_arrow(res.Yref(1:2,k),pos(:,a),':',[0 0.6 0]);
        end
        for a = 1:4
            h(end+1) = plot(x(1,a,k),x(2,a,k),'o','LineWidth',2,'MarkerSize',10,'Color',colors(a+1,:)); %#ok<AGROW>
        end
        h(end+1) = plot(res.Yref(1,k),res.Yref(2,k),'*','LineWidth',1.5,'MarkerSize',8,'Color',[0.85 0.1 0.1]);
        grid on; axis equal;
        xlim([-3 3]); ylim([-3 3]);   % agents of the d = 1.2 m square lie up to 2.6 m from the centre
        xlabel('$x(t)\;(\mathrm{m})$'); ylabel('$y(t)\;(\mathrm{m})$');
        title(sprintf('Formation snapshot at $t = %.2f~\\mathrm{s}$',t(k)));
        legend([h hl hp],{'Reference path','Agent 1','Agent 2','Agent 3','Agent 4','Virtual leader', ...
            'Communication link','Leader information'},'Location','best');
        style(gca);
        export(f,d,names{q});
    end
    e = x-reshape(res.Yref(:,1:T),3,1,T)-res.cfg.offsets;
    comp = {'x','y','\theta'};
    files = {'Fig8_ex_errors','Fig9_ey_errors','Fig11_etheta_errors'};
    titles = {'Formation tracking errors in $x$-direction','Formation tracking errors in $y$-direction', ...
        'Heading-angle tracking errors of all agents vs.\ time'};
    units = {'m','m','rad'}; zoom_units = {'mm','mm','mrad'};
    for c = 1:3
        f = newfig(16,12);
        ec = squeeze(e(c,:,:));
        plot(t,ec.','LineWidth',1.5); grid on;
        legend(arrayfun(@(a) sprintf('$e_{%s_%d}(t)$',comp{c},a),1:4,'UniformOutput',false), ...
            'Location','southeast');
        xlabel('Time (s)'); ylabel(sprintf('$e_%s(t)\\;(\\mathrm{%s})$',comp{c},units{c}));
        title(titles{c});
        lo = min(ec(:)); hi = max(ec(:)); span = hi-lo;
        ylim([lo-0.08*span, max(hi,0)+0.12*span]);
        style(gca);
        inset(f,gca,t,1e3*ec,zoom_units{c});
        export(f,d,files{c});
    end
end

if any(strcmp(which,{'montecarlo','all'}))
    S = load(fullfile(data,'montecarlo_full.mat'));
    rows = S.rows(cellfun(@(r) r.ok,S.rows));
    rows = [rows{:}];
    conds = {'nominal','noise','slip_mismatch','delay_loss','link_failure','link_failure_nosensing', ...
             'link_failure_early','link_failure_early_nosensing','combined'};
    labels = {'Nominal','Noise','Slip+mism.','Delay+loss','Fail 8 s','Fail 8 s, no sens.', ...
              'Fail 1 s','Fail 1 s, no sens.','Combined'};
    present = ismember(conds,cellstr(unique([rows.condition])));
    conds = conds(present); labels = labels(present);
    % Error statistics over the converged runs (final error < 2 cm), solve
    % times over all runs.
    rms = []; group = []; solve = []; sgroup = []; counted = labels;
    for c = 1:numel(conds)
        sel = rows(strcmp([rows.condition],conds{c}));
        conv = [sel.final_position_error] < 0.02;
        rms = [rms, 1e3*[sel(conv).rms_position_error_after]]; %#ok<AGROW>
        group = [group, c*ones(1,nnz(conv))]; %#ok<AGROW>
        solve = [solve, [sel.solve_mean_ms]]; %#ok<AGROW>
        sgroup = [sgroup, c*ones(1,numel(sel))]; %#ok<AGROW>
        counted{c} = sprintf('%s (%d)',labels{c},nnz(conv));
    end
    f = newfig(16,12);
    boxplot(rms,group,'Labels',counted,'Symbol','k.');
    set(gca,'YScale','log','XTickLabelRotation',35); grid on;
    v = rms(rms>0); ylim([10^floor(log10(min(v))) 2*max(v)]);
    ylabel('RMS formation error for $t\ge10$ s (mm)');
    title('RMS formation error of the converged runs');
    style(gca);
    export(f,fullfile(img,'MonteCarlo'),'montecarlo_rms');
    f = newfig(16,12);
    boxplot(solve,sgroup,'Labels',labels,'Symbol','k.');
    set(gca,'YScale','log','XTickLabelRotation',35); hold on; grid on;
    ylim([0.9*min(solve) 1.1*max(solve)]);
    yline(10,'k:','$T_s$','Interpreter','latex','LabelHorizontalAlignment','left','LineWidth',1);
    ylabel('Mean local solve time (ms)');
    title('Mean local solve time of each run');
    style(gca);
    export(f,fullfile(img,'MonteCarlo'),'montecarlo_solve');
end
end

%% ------------------------------------------------------------------------
function f = newfig(w,h)
f = figure('Visible','off','Units','centimeters','Position',[2 2 w h],'Color','w');
end

function style(ax)
set(ax,'FontName','Times New Roman','FontSize',11,'LineWidth',1);
end

function export(f,folder,name)
if ~exist(folder,'dir'), mkdir(folder); end
print(f,fullfile(folder,name),'-depsc','-r600','-painters');
stale = fullfile(folder,[name '-eps-converted-to.pdf']);
if exist(stale,'file'), delete(stale); end
close(f);
end

function inset(f,ax,t,values,unit)
% Magnified view of t >= 10 s placed in the larger empty band of the axes:
% between the legend and the zero line, or between the zero line and the top.
drawnow;
lgd = findobj(f,'Type','Legend');
pos = ax.Position; yl = ylim(ax);
zero_fig = pos(2)+pos(4)*(0-yl(1))/diff(yl);
below = [lgd.Position(2)+lgd.Position(4)+0.04, zero_fig-0.05];
above = [zero_fig+0.06, pos(2)+pos(4)-0.06];
if diff(above) > diff(below), band = above; else, band = below; end
height = min(diff(band),0.34);
box_ = [0.45, band(1)+(diff(band)-height)/2, 0.42, height];
window = t >= 10;
iax = axes(f,'Position',box_,'NextPlot','add','Box','on');
plot(iax,t(window),values(:,window).','LineWidth',1.2);
grid(iax,'on'); xlim(iax,[10 t(end)]);
v = values(:,window); span = max(v(:))-min(v(:)); if span==0, span = 1; end
ylim(iax,[min(v(:))-0.15*span, max(v(:))+0.15*span]);
ylabel(iax,['(' unit ')']);
set(iax,'FontName','Times New Roman','FontSize',9,'LineWidth',0.8);
end

function h = link_arrow(from,to,line_style,color,bend)
% Arrow from position `from` to position `to`, shortened so that it ends
% outside the agent markers, with a filled head of length 0.2 m. For a
% nonzero `bend`, the arrow is a quadratic Bezier curve whose control point
% lies `bend` m to the left of the midpoint.
if nargin < 5, bend = 0; end
v = to-from; len = norm(v);
if len < 0.5, h = plot(nan,nan,line_style,'Color',color,'LineWidth',1.6); return; end
c = (from+to)/2+bend*[-v(2); v(1)]/len;
s = linspace(0,1,401);
P = from*(1-s).^2+c*(2*(1-s).*s)+to*s.^2;
L = [0 cumsum(sqrt(sum(diff(P,1,2).^2,1)))];
body = interp1(L',P',linspace(0.16,L(end)-0.36,100)')';
tip = interp1(L',P',L(end)-0.16)'; base = body(:,end);
u = (tip-base)/norm(tip-base); n = [-u(2); u(1)];
h = plot(body(1,:),body(2,:),line_style,'Color',color,'LineWidth',1.6);
head = [tip, base+0.07*n, base-0.07*n];
patch(head(1,:),head(2,:),color,'EdgeColor',color,'HandleVisibility','off');
end

function pairs = find_pairs(Ad)
[i,j] = find(Ad~=0);
pairs = [i j];
end
