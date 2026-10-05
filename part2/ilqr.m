function [X,U,K,costHistory,info] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fmax,options)
%ILQR Compute a force-limited cart-pole swing-up trajectory.
% The implementation follows the recursions supplied in the project brief,
% clips every force before nonlinear propagation, uses backtracking line
% search, and returns final tracking gains along the accepted trajectory.

if nargin < 9
    options = struct();
end
options = with_default(options,'maxIterations',120);
options = with_default(options,'relativeCostTolerance',1e-5);
options = with_default(options,'feedforwardTolerance',1e-4);
options = with_default(options,'verbose',true);
options = with_default(options,'lineSearch', ...
    [1,0.5,0.25,0.1,0.05,0.01,0.005,0.001]);
options = with_default(options,'minimumRegularization',1e-8);
options = with_default(options,'maximumRegularization',1e8);
options = with_default(options,'maximumAttempts',options.maxIterations);
options = with_default(options,'convergenceRepeats',3);

x0 = x0(:);
U0 = min(Fmax,max(-Fmax,U0(:).'));
N = numel(U0);
if N < 1
    error('CartPole:EmptyHorizon','The iLQR horizon must contain a sample.');
end

% An alpha of zero and zero gains simply rolls out the initial force guess.
[X,U] = forward_pass(x0,zeros(4,N+1),U0,zeros(1,N), ...
    zeros(1,4,N),0,model,Ts,Fmax);
cost = trajectory_cost(X,U,Q,R,Qf);
if ~isfinite(cost)
    error('CartPole:InitialRollout','The initial rollout has nonfinite cost.');
end

costHistory = cost;
acceptedAlphas = zeros(1,options.maxIterations);
regularization = 0;
status = "maximum iterations reached";
acceptedIterations = 0;
attempt = 0;
maximumAttempts = options.maximumAttempts;
smallStepCount = 0;

while acceptedIterations < options.maxIterations && attempt < maximumAttempts
    attempt = attempt+1;
    [A,B] = linearize_trajectory(X,U,model,Ts);

    try
        [d,Kcandidate] = backward_lqr(X,U,A,B,Q,R,Qf,regularization);
    catch backwardError
        regularization = increase_regularization(regularization,options);
        if regularization > options.maximumRegularization
            rethrow(backwardError);
        end
        continue
    end

    % Strong regularization can make d artificially small, so only regard
    % the feedforward term as a convergence test on an essentially
    % unregularized backward pass.
    projectedDirection = min(Fmax,max(-Fmax,U+d))-U;
    if regularization <= options.minimumRegularization ...
            && max(abs(projectedDirection)) < options.feedforwardTolerance
        status = "feedforward tolerance reached";
        break
    end

    accepted = false;
    for alpha = options.lineSearch
        [candidateX,candidateU] = forward_pass(x0,X,U,d,Kcandidate, ...
            alpha,model,Ts,Fmax);
        candidateCost = trajectory_cost(candidateX,candidateU,Q,R,Qf);
        if isfinite(candidateCost) && candidateCost < cost
            accepted = true;
            break
        end
    end

    if ~accepted
        regularization = increase_regularization(regularization,options);
        if regularization > options.maximumRegularization
            status = "line search failed";
            break
        end
        continue
    end

    oldCost = cost;
    controlStep = max(abs(candidateU-U));
    acceptedRegularization = regularization;
    X = candidateX;
    U = candidateU;
    cost = candidateCost;
    costHistory(end+1) = cost; %#ok<AGROW>
    acceptedIterations = acceptedIterations + 1;
    acceptedAlphas(acceptedIterations) = alpha;
    if regularization <= options.minimumRegularization
        regularization = 0;
    else
        regularization = regularization/10;
    end

    relativeDecrease = (oldCost-cost)/max(1,abs(oldCost));
    if options.verbose
        fprintf(['iLQR %3d: J=%11.4f, relative decrease=%8.2e, ' ...
            'max |du|=%8.2e, alpha=%g, lambda=%8.2e\n'], ...
            acceptedIterations,cost,relativeDecrease,controlStep,alpha, ...
            acceptedRegularization);
    end

    % A single small cost change is not convergence on this nonconvex
    % problem: the swing-up routinely crosses shallow plateaus. Require a
    % small control update on several consecutive unregularized steps.
    if options.relativeCostTolerance > 0 ...
            && relativeDecrease < options.relativeCostTolerance ...
            && controlStep < options.feedforwardTolerance ...
            && acceptedRegularization <= options.minimumRegularization
        smallStepCount = smallStepCount+1;
    else
        smallStepCount = 0;
    end
    if smallStepCount >= options.convergenceRepeats
        status = "relative cost tolerance reached";
        break
    end
end

% The project asks for gains recomputed along the final accepted nominal.
[A,B] = linearize_trajectory(X,U,model,Ts);
try
    [~,K] = backward_lqr(X,U,A,B,Q,R,Qf,0);
catch
    [~,K] = backward_lqr(X,U,A,B,Q,R,Qf,max(regularization, ...
        options.minimumRegularization));
end

if attempt >= maximumAttempts && status == "maximum iterations reached"
    status = "attempt limit reached";
end

info.status = char(status);
info.iterations = acceptedIterations;
info.attempts = attempt;
info.finalCost = cost;
info.regularization = regularization;
info.acceptedAlphas = acceptedAlphas(1:acceptedIterations);
info.forceLimit = Fmax;
end

function options = with_default(options,name,value)
if ~isfield(options,name) || isempty(options.(name))
    options.(name) = value;
end
end

function value = increase_regularization(value,options)
if value == 0
    value = options.minimumRegularization;
else
    value = 10*value;
end
end
