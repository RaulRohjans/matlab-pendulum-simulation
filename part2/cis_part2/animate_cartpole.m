%% Animate the cart-pole from the Simscape results (no Mechanics Explorer needed)
% Run AFTER setup_cartpole.m (needs: out, param, Ts, xmax).
% Set saveVideo = true to also write cartpole_swingup.mp4 (or a GIF).

speed     = 0.5;        % 1 = real time, 0.5 = slow motion
tShow     = 5;        % seconds of simulation to animate
saveVideo = false;    % true -> write a video file

% --- data from the simulation
ts = out.x.Time;
xs = squeeze(out.x.Data);
if size(xs,1) == 4 && size(xs,2) ~= 4, xs = xs.'; end      % (samples x 4)
us = squeeze(out.u.Data);  tu = out.u.Time;
p  = xs(:,1);  th = xs(:,3);
L  = param.L;
if exist('P_lqr','var') && exist('c_cap','var') && isfinite(c_cap)
    thw = atan2(sin(th),cos(th));
    V   = sum(([xs(:,1:2) thw xs(:,4)]*P_lqr).*[xs(:,1:2) thw xs(:,4)],2);
    iC  = find(V <= c_cap,1);  if isempty(iC), tCap = inf; else, tCap = ts(iC); end
else
    tCap = 0;
end

% --- figure
fig = figure('Name','Cart-pole animation','Color','w');
ax  = axes(fig); hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
xl  = max(1, max(abs(p))+0.7);
xlim(ax,[-xl xl]); ylim(ax,[-0.7 0.7]);
plot(ax,[-xl+0.1 xl-0.1],[-0.05 -0.05],'Color',[0.6 0.6 0.6],'LineWidth',2);   % track
plot(ax,[-xmax -xmax],[-0.12 0.02],'r','LineWidth',1.5);                       % cart limits
plot(ax,[ xmax  xmax],[-0.12 0.02],'r','LineWidth',1.5);
cart = rectangle(ax,'Position',[p(1)-0.1 -0.05 0.2 0.1], ...
                 'FaceColor',[0.16 0.47 0.84],'EdgeColor','none');
rod  = plot(ax,[0 0],[0 0],'Color',[0.3 0.3 0.3],'LineWidth',4);
bob  = plot(ax,0,0,'o','MarkerSize',10,'MarkerFaceColor',[0.92 0.41 0.2], ...
            'MarkerEdgeColor','none');
txt  = text(ax,-xl+0.05,0.62,'','FontSize',11,'VerticalAlignment','top');

if saveVideo
    vw = VideoWriter('cartpole_swingup','MPEG-4'); vw.FrameRate = 30; open(vw);
end

% --- animation loop (30 frames per second of simulated time)
tFrames = 0:1/30:min(tShow, ts(end));
tic
for tf = tFrames
    k  = find(ts <= tf, 1, 'last');
    ku = find(tu <= tf, 1, 'last');
    % theta = 0 upright, positive counter-clockwise -> tip x = p - L sin(theta)
    px = p(k);  py = 0.05;
    tipx = px - L*sin(th(k));      tipy = py + L*cos(th(k));
    set(cart,'Position',[px-0.1 -0.05 0.2 0.1]);
    set(rod,'XData',[px tipx],'YData',[py tipy]);
    set(bob,'XData',px - L/2*sin(th(k)),'YData',py + L/2*cos(th(k)));   % point mass at L/2
    if tf < tCap, mode = 'swing-up (iLQR tracking)'; else, mode = 'upright LQR'; end
    set(txt,'String',sprintf('t = %4.2f s   %s\nu = %+5.1f N', tf, mode, us(ku)));
    drawnow;
    if saveVideo
        writeVideo(vw, getframe(fig));
    else
        pause(max(0, tf/speed - toc));       % keep (scaled) real time
    end
end
if saveVideo, close(vw); fprintf('Saved cartpole_swingup.mp4\n'); end
