function J = trajectory_cost(X,U,Q,R,Qf)
%TRAJECTORY_COST Quadratic finite-horizon objective from the project brief.
%   J = 1/2*xN'*Qf*xN
%       + 1/2*sum(xk'*Q*xk + uk'*R*uk),  k = 0,...,N-1.

U = U(:).';
N = numel(U);
assert(isequal(size(X), [4, N+1]), ...
    'X must have size 4-by-(numel(U)+1).');

J = 0.5*(X(:,end)'*Qf*X(:,end));
for k = 1:N
    J = J + 0.5*(X(:,k)'*Q*X(:,k) + U(k)'*R*U(k));
end

J = real(J);
end
