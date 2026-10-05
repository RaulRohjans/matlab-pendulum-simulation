function J = trajectory_cost(X,U,Q,R,Qf)
% J is the total cost corresponding to the control trajectory U from
% initial state x0
% Use the objective function given in the project manual
%
%   J = 1/2 x[N]' Qf x[N] + 1/2 sum_{k=0}^{N-1} ( x[k]' Q x[k] + u[k]' R u[k] )
%
% X is 4 x (N+1) (column k+1 holds x[k]), U is 1 x N.

N = numel(U);

%%%% Implement cost calculation here %%%%
xN = X(:,N+1);
J  = 0.5*(xN'*Qf*xN);                 % terminal cost

for k = 1:N                           % running cost, MATLAB index k <-> time k-1
    xk = X(:,k);
    uk = U(k);
    J  = J + 0.5*(xk'*Q*xk + uk'*R*uk);
end

end
