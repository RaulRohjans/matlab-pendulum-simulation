function [Xnew,Unew] = forward_pass(x0,X,U,d,K,alpha,model,Ts,Fmax)
%FORWARD_PASS Roll out one candidate iLQR trajectory on the nonlinear model.
% The project control law is
%   unew = unom + alpha*d - K*(xnew-xnominal).
% Force saturation is applied before each nonlinear propagation.

U = U(:).';
d = d(:).';
N = numel(U);
assert(isequal(size(X), [4, N+1]), 'Unexpected nominal trajectory size.');
assert(numel(d) == N, 'd must contain one entry per control sample.');
assert(isequal(size(K), [1, 4, N]), 'Unexpected feedback gain size.');

Xnew = zeros(4,N+1);
Unew = zeros(1,N);
Xnew(:,1) = x0(:);

for k = 1:N
    stateError = Xnew(:,k) - X(:,k);
    candidateForce = U(k) + alpha*d(k) - K(:,:,k)*stateError;
    Unew(k) = min(Fmax,max(-Fmax,candidateForce));
    Xnew(:,k+1) = discrete_step(Xnew(:,k),Unew(k),model,Ts);
end

end
