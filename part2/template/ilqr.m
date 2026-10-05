function [X,U,K,costHistory] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fmax)

x0 = x0(:);
U0 = min(Fmax,max(-Fmax,U0(:)));
N = numel(U0);

% Initial trajectory
[X,U] = forward_pass(x0,zeros(4,N+1),U0,zeros(1,N),zeros(1,4,N), 0,model,Ts,Fmax);

%%% Implement iLQR here %%%
%%%%% Make sure the implementation has the following:
%%%%% Line search (iterate over different values of alpha in the forward pass) 
%%%%% Stopping and convergence criteria

end