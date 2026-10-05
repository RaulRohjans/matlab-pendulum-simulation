function [X,U,K,costHistory,info] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fmax,opts)
% iLQR swing-up optimization.
%
% Optional opts struct fields (defaults in brackets):
%   maxIter  [200]   maximum number of accepted iterations
%   tol      [1e-4]  stop when relative cost decrease < tol
%   boxed    [true]  use the input-limit-aware backward pass (see backward_lqr)
%   verbose  [false] print progress
%
% Extra output info: .status ('converged' | 'max_iter' | 'line_search_failed'),
%                    .iterations, .mu (final regularization)

if nargin < 9, opts = struct(); end
if ~isfield(opts,'maxIter'), opts.maxIter = 200;   end
if ~isfield(opts,'tol'),     opts.tol     = 1e-4;  end
if ~isfield(opts,'boxed'),   opts.boxed   = true;  end
if ~isfield(opts,'verbose'), opts.verbose = false; end

x0 = x0(:);
U0 = min(Fmax,max(-Fmax,U0(:)));
N = numel(U0);

% Initial trajectory
[X,U] = forward_pass(x0,zeros(4,N+1),U0,zeros(1,N),zeros(1,4,N), 0,model,Ts,Fmax);

%%% Implement iLQR here %%%
%%%%% Make sure the implementation has the following:
%%%%% Line search (iterate over different values of alpha in the forward pass)
%%%%% Stopping and convergence criteria

if opts.boxed, Fbox = Fmax; else, Fbox = inf; end

J = trajectory_cost(X,U,Q,R,Qf);
costHistory = J;

alphas = 0.5.^(0:10);          % line search: 1, 1/2, ..., 1/1024
mu     = 0;                    % regularization (Levenberg-Marquardt)
muMin  = 1e-6;  muMax = 1e8;
status = 'max_iter';
iter   = 0;
relinearize = true;

while iter < opts.maxIter
    % 1) linearize the nonlinear dynamics along the current nominal
    if relinearize
        [A,B] = linearize_trajectory(X,U,model,Ts);
    end

    % 2) backward pass: gains K_k and feedforward d_k
    [d,Kfb] = backward_lqr(X,U,A,B,Q,R,Qf,mu,Fbox);

    % 3) forward pass with line search on alpha
    accepted = false;
    for a = alphas
        [Xn,Un] = forward_pass(x0,X,U,d,Kfb,a,model,Ts,Fmax);
        Jn = trajectory_cost(Xn,Un,Q,R,Qf);
        if Jn < J
            accepted = true;
            break
        end
    end

    if ~accepted
        % no step size decreased the cost: increase regularization and
        % retry the backward pass (same linearization). Give up when mu
        % hits its ceiling -> we are at a (local) minimum / stuck.
        if mu >= muMax
            status = 'line_search_failed';
            break
        end
        mu = max(muMin, 10*mu);
        relinearize = false;
        continue
    end

    % accept the step
    iter = iter + 1;
    relDecrease = (J - Jn)/abs(J);
    X = Xn;  U = Un;  J = Jn;
    costHistory(end+1) = J; %#ok<AGROW>
    relinearize = true;
    if mu > muMin, mu = mu/10; else, mu = 0; end

    if opts.verbose
        fprintf('iLQR it %3d  J = %10.3f  alpha = %6.4f  mu = %.1e\n', iter, J, a, mu);
    end

    % 4) convergence check
    if relDecrease < opts.tol
        status = 'converged';
        break
    end
end

% 5) final feedback gains along the accepted nominal trajectory
%    (unconstrained recursion, no regularization) used for tracking
[A,B] = linearize_trajectory(X,U,model,Ts);
[~,K] = backward_lqr(X,U,A,B,Q,R,Qf);

info.status     = status;
info.iterations = iter;
info.mu         = mu;
if opts.verbose
    fprintf('iLQR finished: %s after %d iterations, J = %.3f\n', status, iter, J);
end

end
