function study = capture_zone_experiment(makePlots)
%CAPTURE_ZONE_EXPERIMENT Empirically determine a saturated-LQR capture zone.
% The nonlinear plant is sampled in all four state dimensions. The proposed
% criterion is the conservative Riccati level set x'*P*x <= c_cap, where
% c_cap is 75 percent of the smallest sampled failure level.

if nargin < 1
    makePlots = true;
end

projectRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(projectRoot,'results');
if makePlots && ~isfolder(resultsDir)
    mkdir(resultsDir);
end

model = struct('M',1.0,'m',1.0,'ell',0.25,'bc',0.1,'bp',0.05,'g',9.81);
Ts = 0.01;
Fmax = 20;
xmax = 0.5;
[A,B] = upright_linear_model(model);
transition = expm([A,B;zeros(1,5)]*Ts);
Ad = transition(1:4,1:4);
Bd = transition(1:4,5);
Q_lqr = diag([1/0.5^2,1/1.0^2,1/0.2^2,1/2.0^2]);
R_lqr = 1/Fmax^2;
[K_lqr,P_lqr] = discrete_lqr(Ad,Bd,Q_lqr,R_lqr);

% Independent signed sweeps establish readable one-state-at-a-time limits.
axisRanges = {0:0.01:0.5,0:0.05:3,0:0.01:1.3,0:0.05:7};
axisNames = {'position_m','velocity_m_s','angle_rad','angular_velocity_rad_s'};
axisLimits = zeros(1,4);
for stateIndex = 1:4
    values = axisRanges{stateIndex};
    positive = zeros(4,numel(values));
    positive(stateIndex,:) = values;
    negative = -positive;
    okPositive = simulate_saturated_lqr(positive,K_lqr,model,Ts,Fmax,xmax,6);
    okNegative = simulate_saturated_lqr(negative,K_lqr,model,Ts,Fmax,xmax,6);
    bothDirections = okPositive & okNegative;
    firstFailure = find(~bothDirections,1);
    if isempty(firstFailure)
        axisLimits(stateIndex) = values(end);
    else
        axisLimits(stateIndex) = values(max(1,firstFailure-1));
    end
end

% Combined-state Monte Carlo is essential: independent axis limits alone do
% not define a safe four-dimensional capture region.
rng(24,'twister');
nSamples = 6000;
sampledStates = [0.45*(2*rand(1,nSamples)-1);
                 2.5*(2*rand(1,nSamples)-1);
                 1.2*(2*rand(1,nSamples)-1);
                 6.0*(2*rand(1,nSamples)-1)];
sampledSuccess = simulate_saturated_lqr(sampledStates,K_lqr,model, ...
    Ts,Fmax,xmax,6);
level = sum(sampledStates.*(P_lqr*sampledStates),1);
failureLevels = level(~sampledSuccess);
if isempty(failureLevels)
    firstFailureLevel = max(level);
else
    firstFailureLevel = min(failureLevels);
end
c_cap = 0.75*firstFailureLevel;
inside = level <= c_cap;
if ~any(inside) || ~all(sampledSuccess(inside))
    error('CartPole:CaptureStudy','Could not construct a sampled-safe level set.');
end

study.c_cap = c_cap;
study.firstFailureLevel = firstFailureLevel;
study.axisLimits = axisLimits;
study.axisNames = axisNames;
study.sampleCount = nSamples;
study.sampledSuccessRate = mean(sampledSuccess);
study.insideCount = nnz(inside);
study.insideSuccessRate = mean(sampledSuccess(inside));
study.K_lqr = K_lqr;
study.P_lqr = P_lqr;

fprintf('\nEmpirical saturated-LQR capture study\n');
fprintf('  sampled states: %d, overall stabilized: %.1f%%\n', ...
    nSamples,100*study.sampledSuccessRate);
fprintf('  proposed criterion: x''P x <= %.3f\n',c_cap);
fprintf('  validation inside criterion: %d/%d stabilized\n', ...
    nnz(sampledSuccess(inside)),study.insideCount);
for stateIndex = 1:4
    fprintf('  signed-axis limit %-24s %.3f\n', ...
        axisNames{stateIndex},axisLimits(stateIndex));
end

if makePlots
    captureSummary = table(c_cap,firstFailureLevel,nSamples, ...
        study.sampledSuccessRate,study.insideCount,study.insideSuccessRate, ...
        axisLimits(1),axisLimits(2),axisLimits(3),axisLimits(4), ...
        'VariableNames',{'CaptureLevel','FirstFailureLevel','SampleCount', ...
        'OverallSuccessRate','InsideCount','InsideSuccessRate', ...
        'PositionLimit_m','VelocityLimit_m_s','AngleLimit_rad', ...
        'AngularVelocityLimit_rad_s'});
    writetable(captureSummary,fullfile(resultsDir,'capture_zone_summary.csv'));

    thetaValues = linspace(-1.2,1.2,71);
    omegaValues = linspace(-6,6,71);
    [thetaGrid,omegaGrid] = meshgrid(thetaValues,omegaValues);
    sliceStates = [zeros(2,numel(thetaGrid));thetaGrid(:)';omegaGrid(:)'];
    thetaOmegaSuccess = reshape(simulate_saturated_lqr(sliceStates,K_lqr, ...
        model,Ts,Fmax,xmax,6),size(thetaGrid));

    positionValues = linspace(-0.5,0.5,71);
    velocityValues = linspace(-2.5,2.5,71);
    [positionGrid,velocityGrid] = meshgrid(positionValues,velocityValues);
    sliceStates = [positionGrid(:)';velocityGrid(:)';zeros(2,numel(positionGrid))];
    positionVelocitySuccess = reshape(simulate_saturated_lqr(sliceStates, ...
        K_lqr,model,Ts,Fmax,xmax,6),size(positionGrid));

    fig = figure('Name','Part 2 - empirical LQR capture zone','Color','w');
    tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    firstAxes = nexttile;
    imagesc(firstAxes,thetaValues,omegaValues,double(thetaOmegaSuccess));
    axis(firstAxes,'xy'); hold(firstAxes,'on'); grid(firstAxes,'on');
    clim(firstAxes,[0,1]);
    colormap(firstAxes,[0.88,0.88,0.88;0.35,0.68,0.42]);
    contour(thetaGrid,omegaGrid,quadratic_level(thetaGrid,omegaGrid, ...
        P_lqr([3,4],[3,4])),[c_cap,c_cap],'k','LineWidth',1.5);
    xlabel('\theta [rad]'); ylabel('\omega [rad/s]');
    title('\theta-\omega slice (p=v=0)');
    secondAxes = nexttile;
    imagesc(secondAxes,positionValues,velocityValues, ...
        double(positionVelocitySuccess));
    axis(secondAxes,'xy'); hold(secondAxes,'on'); grid(secondAxes,'on');
    clim(secondAxes,[0,1]);
    colormap(secondAxes,[0.88,0.88,0.88;0.35,0.68,0.42]);
    contour(positionGrid,velocityGrid,quadratic_level(positionGrid,velocityGrid, ...
        P_lqr([1,2],[1,2])),[c_cap,c_cap],'k','LineWidth',1.5);
    xlabel('p [m]'); ylabel('v [m/s]');
    title('p-v slice (\theta=\omega=0)');
    exportgraphics(fig,fullfile(resultsDir,'01_capture_zone.png'),'Resolution',160);
end
end

function values = quadratic_level(a,b,P)
values = P(1,1)*a.^2+2*P(1,2)*a.*b+P(2,2)*b.^2;
end

function [A,B] = upright_linear_model(p)
J = p.m*p.ell^2;
Delta = (p.M+p.m)*J-(p.m*p.ell)^2;
A = [0,1,0,0;
     0,-J*p.bc/Delta,p.m^2*p.g*p.ell^2/Delta,-p.m*p.ell*p.bp/Delta;
     0,0,0,1;
     0,-p.m*p.ell*p.bc/Delta,(p.M+p.m)*p.m*p.g*p.ell/Delta, ...
        -(p.M+p.m)*p.bp/Delta];
B = [0;J/Delta;0;p.m*p.ell/Delta];
end

function [K,P] = discrete_lqr(A,B,Q,R)
if ~isempty(which('dlqr'))
    [K,P] = dlqr(A,B,Q,R);
    return
end
P = Q;
for iteration = 1:100000
    K = (R+B'*P*B)\(B'*P*A);
    nextP = Q+A'*P*(A-B*K);
    if norm(nextP-P,'fro') <= 1e-11*max(1,norm(P,'fro'))
        P = nextP;
        break
    end
    P = nextP;
end
K = (R+B'*P*B)\(B'*P*A);
end
