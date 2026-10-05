function [experimentTable,workflowResult] = run_part2_experiments(makePlots)
%RUN_PART2_EXPERIMENTS Required iLQR sensitivity and feasibility studies.
% Covers force/cart limits, cost weights, horizon, and two initial guesses.
% Clicking Run (or calling this function without inputs) first runs the
% complete Part 2 workflow, so all four result figures remain visible:
% capture zone, nominal iLQR, Simscape validation, and experiments.
% Use run_part2_experiments(false) for experiments without figures or
% Simulink.

if nargin < 1
    makePlots = true;
end
projectRoot = fileparts(mfilename('fullpath'));
resultsDir = fullfile(projectRoot,'results');
if ~isfolder(resultsDir),mkdir(resultsDir);end

workflowResult = struct();
if makePlots
    % Keep figures visible in ordinary MATLAB windows when this file is run
    % from the Editor. setup_cartpole also exports the same figures to the
    % results directory.
    set(groot,'defaultFigureVisible','on');
    if usejava('desktop')
        set(groot,'defaultFigureWindowStyle','normal');
    end
    workflowResult = setup_cartpole(struct( ...
        'runSimulink',true, ...
        'makePlots',true, ...
        'runCaptureStudy',true, ...
        'showMultibodyExplorer',usejava('desktop')));
    capture = workflowResult.captureStudy;
else
    capture = capture_zone_experiment(false);
end

model = struct('M',1.0,'m',1.0,'ell',0.25,'bc',0.1,'bp',0.05,'g',9.81);
Ts = 0.01;
x0 = [0;0;pi;0];
P_lqr = capture.P_lqr;
c_cap = capture.c_cap;

base = struct('Fmax',20,'xmax',0.5,'T',2.0,'Qp',10, ...
    'R',1e-2,'QfScale',1,'guess','zero');
scenarios = repmat(base,10,1);
labels = ["baseline";"force 10 N";"force 30 N";"cart limit 0.30 m"; ...
    "cart limit 0.70 m";"horizon 1.5 s";"horizon 3.0 s"; ...
    "larger cart weight";"smaller force weight";"sinusoidal initial guess"];
categories = ["baseline";"force";"force";"cart limit";"cart limit"; ...
    "horizon";"horizon";"weights";"weights";"initial guess"];
scenarios(2).Fmax = 10;
scenarios(3).Fmax = 30;
scenarios(4).xmax = 0.30;
scenarios(5).xmax = 0.70;
scenarios(6).T = 1.5;
scenarios(7).T = 3.0;
scenarios(8).Qp = 50;
scenarios(9).R = 1e-3;
scenarios(10).guess = 'sine';

n = numel(scenarios);
accepted = false(n,1);
iterations = zeros(n,1);
finalCost = nan(n,1);
terminalLevel = nan(n,1);
maximumPosition = nan(n,1);
maximumForce = nan(n,1);
status = strings(n,1);

for trial = 1:n
    s = scenarios(trial);
    N = round(s.T/Ts);
    time = (0:N-1)*Ts;
    if strcmp(s.guess,'sine')
        U0 = 0.35*s.Fmax*sin(2*pi*time);
    else
        U0 = zeros(1,N);
    end
    Q = diag([s.Qp,0.1,0.1,0.1]);
    Qf = s.QfScale*1000*diag([1,1,10,1]);
    opts = struct('maxIterations',120,'relativeCostTolerance',1e-6, ...
        'feedforwardTolerance',1e-4,'verbose',false);
    try
        [X,U,~,costs,info] = ilqr(x0,U0,model,Ts,Q,s.R,Qf, ...
            0.9*s.Fmax,opts);
        terminal = X(:,end);
        terminal(3) = atan2(sin(terminal(3)),cos(terminal(3)));
        terminalLevel(trial) = terminal'*P_lqr*terminal;
        maximumPosition(trial) = max(abs(X(1,:)));
        maximumForce(trial) = max(abs(U));
        finalCost(trial) = costs(end);
        iterations(trial) = info.iterations;
        status(trial) = string(info.status);
        accepted(trial) = terminalLevel(trial) <= c_cap ...
            && maximumPosition(trial) <= s.xmax+1e-9 ...
            && maximumForce(trial) <= 0.9*s.Fmax+1e-9 ...
            && all(diff(costs) < 1e-10);
    catch trialError
        status(trial) = "error: "+string(trialError.identifier);
    end
    fprintf('%-27s %-6s  J=%9.2f  V_N=%8.2f  max|p|=%.3f  max|u|=%.2f\n', ...
        labels(trial),pass_fail(accepted(trial)),finalCost(trial), ...
        terminalLevel(trial),maximumPosition(trial),maximumForce(trial));
end

Fmax = reshape([scenarios.Fmax],[],1);
xmax = reshape([scenarios.xmax],[],1);
horizon = reshape([scenarios.T],[],1);
cartWeight = reshape([scenarios.Qp],[],1);
forceWeight = reshape([scenarios.R],[],1);
initialGuess = strings(n,1);
for trial = 1:n
    initialGuess(trial) = scenarios(trial).guess;
end
experimentTable = table(categories,labels,Fmax,xmax,horizon, ...
    cartWeight,forceWeight,initialGuess,accepted,status,iterations,finalCost, ...
    terminalLevel,maximumPosition,maximumForce,'VariableNames', ...
    {'Category','Scenario','Fmax_N','xmax_m','Horizon_s','CartWeight', ...
     'ForceWeight','InitialGuess','Accepted','Status','Iterations', ...
     'FinalCost','TerminalCaptureLevel','MaxPosition_m','MaxForce_N'});
writetable(experimentTable,fullfile(resultsDir,'ilqr_experiments.csv'));

if makePlots
    fig = figure('Name','Part 2 - iLQR experiments','Color','w');
    tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
    nexttile;
    bar(categorical(experimentTable.Scenario),experimentTable.FinalCost);
    ylabel('Final cost'); grid on; title('iLQR sensitivity experiments');
    nexttile;
    bar(categorical(experimentTable.Scenario),experimentTable.MaxPosition_m);
    ylabel('Maximum |p| [m]'); grid on;
    exportgraphics(fig,fullfile(resultsDir,'04_ilqr_experiments.png'), ...
        'Resolution',160);
    drawnow;
end
end

function value = pass_fail(condition)
if condition,value='ACCEPT';else,value='REJECT';end
end
