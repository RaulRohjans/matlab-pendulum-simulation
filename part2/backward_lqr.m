function [d,K] = backward_lqr(X,U,A,B,Q,R,Qf,regularization)
%BACKWARD_LQR Backward recursion for the LQ tracking subproblem.
% d is the feedforward term called k in the project brief. The different
% name avoids confusion with MATLAB's time-step index. The optional scalar
% regularization is added to R + B'*S*B when iLQR needs a safer step.

if nargin < 8 || isempty(regularization)
    regularization = 0;
end

U = U(:).';
N = numel(U);
assert(isequal(size(X), [4, N+1]), 'Unexpected state trajectory size.');
assert(isequal(size(A), [4, 4, N]), 'Unexpected A trajectory size.');
assert(isequal(size(B), [4, 1, N]), 'Unexpected B trajectory size.');

d = zeros(1,N);
K = zeros(1,4,N);

% Boundary conditions in Section 4 of the brief.
S = Qf;
s = Qf*X(:,N+1);

for k = N:-1:1
    Ak = A(:,:,k);
    Bk = B(:,:,k);
    H = R + Bk'*S*Bk + regularization;

    if ~isfinite(H) || H <= eps
        error('CartPole:BackwardPass', ...
            'The control Hessian is not positive at sample %d.', k);
    end

    Kk = H \ (Bk'*S*Ak);
    dk = -H \ (R*U(k) + Bk'*s);

    % Use the old S and s on both right-hand sides before stepping back.
    sNew = Q*X(:,k) + Ak'*s + Ak'*S*Bk*dk;
    SNew = Q + Ak'*S*(Ak - Bk*Kk);

    d(k) = dk;
    K(:,:,k) = Kk;
    s = sNew;
    S = 0.5*(SNew + SNew');
end
end
