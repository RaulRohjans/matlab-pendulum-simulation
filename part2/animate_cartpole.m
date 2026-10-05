function animate_cartpole(result,speed)
%ANIMATE_CARTPOLE Display the Simscape result without Multibody Explorer.
% Usage after setup_cartpole:
%   result = setup_cartpole;
%   animate_cartpole(result)

if nargin < 2 || isempty(speed),speed=1;end
if ~isfield(result,'simulation') || ~result.simulation.ran
    error('CartPole:NoSimulation','Run setup_cartpole with Simulink enabled first.');
end

output = result.simulation.output;
stateSeries = output.get('x');
forceSeries = output.get('u');
time = stateSeries.Time(:);
state = squeeze(stateSeries.Data);
if size(state,2) ~= 4 && size(state,1) == 4,state=state.';end
force = squeeze(forceSeries.Data);force=force(:);
forceTime = forceSeries.Time(:);

poleLength = result.parameters.L;
xmax = result.xmax;
fig = figure('Name','Part 2 - cart-pole animation','Color','w');
ax = axes(fig);hold(ax,'on');axis(ax,'equal');axis(ax,'off');
horizontalLimit = max(1,max(abs(state(:,1)))+0.65);
xlim(ax,[-horizontalLimit,horizontalLimit]);ylim(ax,[-0.65,0.65]);
plot(ax,[-horizontalLimit,horizontalLimit],[-0.05,-0.05], ...
    'Color',[0.55,0.55,0.55],'LineWidth',2);
plot(ax,[-xmax,-xmax],[-0.12,0.03],'r--');
plot(ax,[xmax,xmax],[-0.12,0.03],'r--');
cart = rectangle(ax,'Position',[-0.12,-0.05,0.24,0.10], ...
    'FaceColor',[0.14,0.45,0.82],'EdgeColor','none');
rod = plot(ax,[0,0],[0,poleLength],'Color',[0.2,0.2,0.2],'LineWidth',4);
bob = plot(ax,0,poleLength,'o','MarkerSize',11, ...
    'MarkerFaceColor',[0.92,0.36,0.18],'MarkerEdgeColor','none');
caption = text(ax,-horizontalLimit+0.04,0.60,'','VerticalAlignment','top');

frameTimes = 0:1/30:time(end);
animationStart = tic;
for frameTime = frameTimes
    stateIndex = find(time<=frameTime,1,'last');
    forceIndex = find(forceTime<=frameTime,1,'last');
    position = state(stateIndex,1);
    theta = state(stateIndex,3);
    pivotHeight = 0.05;
    tipX = position-poleLength*sin(theta);
    tipY = pivotHeight+poleLength*cos(theta);
    set(cart,'Position',[position-0.12,-0.05,0.24,0.10]);
    set(rod,'XData',[position,tipX],'YData',[pivotHeight,tipY]);
    set(bob,'XData',tipX,'YData',tipY);
    set(caption,'String',sprintf('t = %.2f s    force = %+.2f N', ...
        frameTime,force(forceIndex)));
    drawnow;
    pause(max(0,frameTime/speed-toc(animationStart)));
end
end
