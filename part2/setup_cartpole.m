function result = setup_cartpole(options)
%SETUP_CARTPOLE Complete Part 2: iLQR swing-up followed by upright LQR.
% Run without inputs for the complete experiment:
%   setup_cartpole
% For a faster equation-only check:
%   setup_cartpole(struct('runSimulink',false,'makePlots',false))

if nargin < 1
    options = struct();
end
options = add_default(options,'runSimulink',true);
options = add_default(options,'makePlots',true);
options = add_default(options,'verbose',true);
options = add_default(options,'runCaptureStudy',true);
options = add_default(options,'showMultibodyExplorer',usejava('desktop'));

projectRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(projectRoot,'results');
if ~isfolder(resultsDir)
    mkdir(resultsDir);
end

%% Physical plant and initial downward equilibrium
param.m = 1.0;
param.M = 1.0;
param.L = 0.5;
param.bc = 0.1;
param.bp = 0.05;
param.p0 = 0;
param.v0 = 0;
param.theta0 = 180;             % Simscape uses degrees; 0 is upright
param.omega0 = 0;

model.M = param.M;
model.m = param.m;
model.ell = param.L/2;
model.bc = param.bc;
model.bp = param.bp;
model.g = 9.81;
x0 = [param.p0;param.v0;deg2rad(param.theta0);param.omega0];

Ts = 0.01;
Fmax = 20;
xmax = 0.5;

%% Upright linearization, exact ZOH discretization, and saturated LQR
[A,B] = upright_linear_model(model);
zoh = expm([A,B;zeros(1,5)]*Ts);
Ad = zoh(1:4,1:4);
Bd = zoh(1:4,5);

% Bryson-style normalized weights. Angle is weighted most strongly because
% upright balance is the primary goal; position reflects the 0.5 m track.
acceptableState = [0.5,1.0,0.2,2.0];
Q_lqr = diag(1./acceptableState.^2);
R_lqr = 1/Fmax^2;
[K_lqr,P_lqr] = discrete_lqr(Ad,Bd,Q_lqr,R_lqr);
closedLoopRadius = max(abs(eig(Ad-Bd*K_lqr)));
if closedLoopRadius >= 1
    error('CartPole:UnstableLQR','The discrete upright LQR is not stable.');
end

% The empirical experiment varies all states and validates the proposed
% Riccati level set on combined-state samples of the nonlinear plant.
if options.runCaptureStudy
    captureStudy = capture_zone_experiment(options.makePlots);
    c_cap = captureStudy.c_cap;
else
    c_cap = 55; % conservative value reproduced by capture_zone_experiment
    captureStudy = struct('c_cap',c_cap);
end

% Independent nonlinear sanity check inside the capture region.
lqrTestState = [0;0;deg2rad(10);0];
lqrStabilized = simulate_saturated_lqr(lqrTestState,K_lqr,model, ...
    Ts,Fmax,xmax,6);
if ~lqrStabilized
    error('CartPole:LQRValidation','The nonlinear upright LQR check failed.');
end

%% iLQR swing-up optimization
T = 2.0;
N = round(T/Ts);
T = N*Ts;
U0 = zeros(1,N);

% Running cost discourages cart travel but stays light on pole motion so the
% optimizer can swing through pi radians. The terminal angle weight is large
% to place the final state inside the LQR capture region.
Q = diag([10,0.1,0.1,0.1]);
% The experiment sweep showed that 1e-2 avoids a saturated local minimum
% reached with 1e-3 while still using the available force effectively.
R = 1e-2;
Qf = 1000*diag([1,1,10,1]);
planningForceLimit = 0.9*Fmax; % preserve feedback authority in Simscape

ilqrOptions.maxIterations = 150;
ilqrOptions.relativeCostTolerance = 1e-6;
ilqrOptions.feedforwardTolerance = 1e-4;
ilqrOptions.verbose = options.verbose;
[Xnom,Unom,Knom,costHistory,optimization] = ilqr(x0,U0,model,Ts, ...
    Q,R,Qf,planningForceLimit,ilqrOptions);

terminalForCapture = Xnom(:,end);
terminalForCapture(3) = atan2(sin(terminalForCapture(3)), ...
    cos(terminalForCapture(3)));
terminalLevel = terminalForCapture'*P_lqr*terminalForCapture;
nominalChecks.finite = all(isfinite(Xnom),'all') && all(isfinite(Unom));
nominalChecks.capture = terminalLevel <= c_cap;
nominalChecks.cart = max(abs(Xnom(1,:))) <= xmax+1e-9;
nominalChecks.force = max(abs(Unom)) <= planningForceLimit+1e-9;
nominalChecks.monotoneCost = all(diff(costHistory) < 1e-10);
nominalAccepted = all(structfun(@(value) logical(value),nominalChecks));

fprintf('\nNominal iLQR trajectory\n');
fprintf('  status: %s (%d accepted iterations)\n', ...
    optimization.status,optimization.iterations);
fprintf('  cost: %.3f -> %.3f\n',costHistory(1),costHistory(end));
fprintf('  terminal x''P x: %.3f (capture threshold %.3f)\n', ...
    terminalLevel,c_cap);
fprintf('  max |p|: %.3f m (limit %.3f), max |u|: %.3f N\n', ...
    max(abs(Xnom(1,:))),xmax,max(abs(Unom)));
fprintf('  feasibility: %s\n',pass_fail(nominalAccepted));

if ~nominalAccepted
    error('CartPole:NominalRejected', ...
        ['The optimized trajectory failed a feasibility check. Increase the ' ...
         'horizon or terminal weight, or retune the cart-position weight.']);
end

if options.makePlots
    plot_nominal(resultsDir,Xnom,Unom,costHistory,Ts,xmax,Fmax);
end

%% Nonlinear Simscape validation and automatic iLQR-to-LQR switching
simulation = struct('ran',false);
if options.runSimulink
    modelFile = fullfile(projectRoot,'cartpole.slx');
    if ~isfile(modelFile)
        error('CartPole:MissingModel','Missing Simulink model: %s',modelFile);
    end

    load_system(modelFile);
    modelName = 'cartpole';
    % Show the physical cart and pendulum for normal desktop MATLAB. The
    % viewer is disabled automatically in headless MATLAB, where numerical
    % validation and the ordinary MATLAB plots still work.
    if options.showMultibodyExplorer
        explorerSetting = 'on';
    else
        explorerSetting = 'off';
    end
    set_param(modelName,'SimMechanicsOpenEditorOnUpdate',explorerSetting);

    simulationEnd = max(8,T+5);
    simulationInput = Simulink.SimulationInput(modelName);
    variableNames = {'param','Ts','Fmax','K_lqr','Knom','Unom','Xnom', ...
        'P_lqr','c_cap'};
    variableValues = {param,Ts,Fmax,K_lqr,Knom,Unom,Xnom,P_lqr,c_cap};
    for variableIndex = 1:numel(variableNames)
        simulationInput = simulationInput.setVariable( ...
            variableNames{variableIndex},variableValues{variableIndex}, ...
            'Workspace',modelName);
    end
    simulationInput = simulationInput.setModelParameter( ...
        'StopTime',num2str(simulationEnd),'ReturnWorkspaceOutputs','on');
    output = sim(simulationInput);

    stateSeries = output.get('x');
    forceSeries = output.get('u');
    time = stateSeries.Time(:);
    state = squeeze(stateSeries.Data);
    if size(state,2) ~= 4 && size(state,1) == 4
        state = state.';
    end
    force = squeeze(forceSeries.Data);
    force = force(:);
    forceTime = forceSeries.Time(:);

    wrappedState = state;
    wrappedState(:,3) = atan2(sin(state(:,3)),cos(state(:,3)));
    captureLevel = sum((wrappedState*P_lqr).*wrappedState,2);
    captureIndex = find(captureLevel <= c_cap,1);
    if isempty(captureIndex)
        captureTime = NaN;
    else
        captureTime = time(captureIndex);
    end

    finalWindow = time >= time(end)-0.5;
    simulationChecks.captured = ~isnan(captureTime) && captureTime <= T+Ts;
    simulationChecks.settled = all(abs(wrappedState(finalWindow,3)) < deg2rad(3)) ...
        && all(abs(wrappedState(finalWindow,4)) < 0.12);
    simulationChecks.cart = max(abs(state(:,1))) <= xmax+0.02;
    simulationChecks.force = max(abs(force)) <= Fmax+1e-8;
    simulationPassed = all(structfun(@(value) logical(value),simulationChecks));

    fprintf('\nNonlinear Simscape validation\n');
    if isnan(captureTime)
        fprintf('  capture: FAILED before the trajectory ended\n');
    else
        fprintf('  capture time: %.3f s (trajectory horizon %.3f s)\n', ...
            captureTime,T);
    end
    fprintf('  final state: p=%+.4f, v=%+.4f, theta=%+.4f, omega=%+.4f\n', ...
        state(end,1),state(end,2),wrappedState(end,3),state(end,4));
    fprintf('  max |p|=%.3f m, max |u|=%.3f N, validation=%s\n', ...
        max(abs(state(:,1))),max(abs(force)),pass_fail(simulationPassed));

    if options.makePlots
        plot_simscape(resultsDir,time,state,wrappedState,forceTime,force, ...
            Xnom,Unom,Ts,captureTime,xmax,Fmax);
    end

    metrics = table(captureTime,state(end,1),state(end,2), ...
        wrappedState(end,3),state(end,4),max(abs(state(:,1))), ...
        max(abs(force)),simulationPassed,'VariableNames', ...
        {'CaptureTime_s','FinalPosition_m','FinalVelocity_m_s', ...
         'FinalAngle_rad','FinalAngularVelocity_rad_s','MaxPosition_m', ...
         'MaxForce_N','Passed'});
    writetable(metrics,fullfile(resultsDir,'part2_metrics.csv'));

    simulation.ran = true;
    simulation.passed = simulationPassed;
    simulation.captureTime = captureTime;
    simulation.checks = simulationChecks;
    simulation.output = output;
    if ~simulationPassed
        warning('CartPole:SimscapeValidation', ...
            'The Simscape run completed but did not satisfy every acceptance check.');
    end
end

result.parameters = param;
result.model = model;
result.Ts = Ts;
result.Fmax = Fmax;
result.xmax = xmax;
result.Ad = Ad;
result.Bd = Bd;
result.K_lqr = K_lqr;
result.P_lqr = P_lqr;
result.closedLoopRadius = closedLoopRadius;
result.captureStudy = captureStudy;
result.c_cap = c_cap;
result.Xnom = Xnom;
result.Unom = Unom;
result.Knom = Knom;
result.costHistory = costHistory;
result.optimization = optimization;
result.nominalTerminalLevel = terminalLevel;
result.nominalChecks = nominalChecks;
result.simulation = simulation;

fprintf('\nPart 2 workflow complete. Results: %s\n',resultsDir);
end

function options = add_default(options,name,value)
if ~isfield(options,name) || isempty(options.(name))
    options.(name) = value;
end
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

function plot_nominal(resultsDir,X,U,costHistory,Ts,xmax,Fmax)
time = (0:size(X,2)-1)*Ts;
fig = figure('Name','Part 2 - optimized iLQR swing-up','Color','w');
tiledlayout(fig,3,2,'Padding','compact','TileSpacing','compact');
labels = {'p [m]','v [m/s]','\theta [rad]','\omega [rad/s]'};
for stateIndex = 1:4
    nexttile;
    plot(time,X(stateIndex,:),'LineWidth',1.25); grid on;
    ylabel(labels{stateIndex}); xlabel('Time [s]');
    if stateIndex == 1
        hold on; yline([xmax,-xmax],'r--');
    end
end
nexttile;
stairs(time(1:end-1),U,'LineWidth',1.25); grid on; hold on;
yline([Fmax,-Fmax],'r--'); ylabel('Force [N]'); xlabel('Time [s]');
nexttile;
semilogy(0:numel(costHistory)-1,costHistory,'o-','LineWidth',1.1);
grid on; xlabel('Accepted iteration'); ylabel('Cost');
exportgraphics(fig,fullfile(resultsDir,'02_ilqr_nominal.png'),'Resolution',160);
end

function plot_simscape(resultsDir,time,state,wrappedState,forceTime,force, ...
    Xnom,Unom,Ts,captureTime,xmax,Fmax)
nominalTime = (0:size(Xnom,2)-1)*Ts;
fig = figure('Name','Part 2 - Simscape swing-up and stabilization','Color','w');
tiledlayout(fig,3,2,'Padding','compact','TileSpacing','compact');
labels = {'p [m]','v [m/s]','\theta [rad]','\omega [rad/s]'};
for stateIndex = 1:4
    nexttile;
    if stateIndex == 3
        plot(time,unwrap(state(:,3)),'LineWidth',1.2);
    else
        plot(time,state(:,stateIndex),'LineWidth',1.2);
    end
    hold on; grid on;
    plot(nominalTime,Xnom(stateIndex,:),'--','LineWidth',1.0);
    if ~isnan(captureTime),xline(captureTime,'k:');end
    if stateIndex == 1,yline([xmax,-xmax],'r--');end
    ylabel(labels{stateIndex}); xlabel('Time [s]');
end
nexttile;
stairs(forceTime,force,'LineWidth',1.2); hold on; grid on;
stairs(nominalTime(1:end-1),Unom,'--','LineWidth',1.0);
yline([Fmax,-Fmax],'r--');
if ~isnan(captureTime),xline(captureTime,'k:');end
ylabel('Force [N]'); xlabel('Time [s]');
nexttile;
plot(time,wrappedState(:,3),'LineWidth',1.2); hold on; grid on;
yline([deg2rad(3),-deg2rad(3)],'r--');
if ~isnan(captureTime),xline(captureTime,'k:');end
ylabel('Wrapped angle [rad]'); xlabel('Time [s]');
legend('Simscape','Nominal/capture bound','Location','best');
exportgraphics(fig,fullfile(resultsDir,'03_simscape_validation.png'),'Resolution',160);
end

function value = pass_fail(condition)
if condition
    value = 'PASS';
else
    value = 'FAIL';
end
end
