function [d,K] = backward_lqr(X,U,A,B,Q,R,Qf,mu,Fmax)
% A is of size (4,4,N), A matrices corresponding to linearized trajectory
% B is of size (4,1,N), B matrices correspondng to linearized trajectory
%
% Backward recursion of the LQ tracking problem solved in each iLQR
% iteration (project manual, Sec. 4):
%   K_k = (R + B_k' S_{k+1} B_k)^-1 B_k' S_{k+1} A_k
%   k_k = -(R + B_k' S_{k+1} B_k)^-1 (R ubar[k] + B_k' s_{k+1})
%   S_k = Q + A_k' S_{k+1} (A_k - B_k K_k)
%   s_k = Q xbar[k] + A_k' s_{k+1} + A_k' S_{k+1} B_k k_k
% with S_N = Qf, s_N = Qf xbar[N].
%
% Optional inputs (not in the manual, used to make iLQR robust):
%   mu   - Levenberg-Marquardt regularization added to (R + B'SB).
%          Large mu -> small, gradient-like steps. Default 0.
%   Fmax - force limit. If ubar[k] + k_k would violate the limit, k_k is
%          clamped to the limit and K_k is set to zero (the input is
%          saturated, so it cannot react to state deviations). For a scalar
%          input this is the exact solution of the box-constrained QP
%          (control-limited DDP). Default inf (= plain recursion above).

if nargin < 8 || isempty(mu),   mu   = 0;   end
if nargin < 9 || isempty(Fmax), Fmax = inf; end

N = numel(U);  %number of time steps

%%% Initialization
d = zeros(1,N); % feedforward terms
% Note: d corresponds to the feedforward k term in the class notes

K = zeros(1,4,N); %state feedback gains

%%% Implement backward pass here
S = Qf;                 % S_N
s = Qf*X(:,N+1);        % s_N

for k = N:-1:1
    Ak = A(:,:,k);
    Bk = B(:,:,k);

    Quu = R + Bk'*S*Bk + mu;            % (R + B'SB), scalar here
    Kk  = Quu \ (Bk'*S*Ak);             % 1x4 feedback gain
    dk  = -Quu \ (R*U(k) + Bk'*s);      % feedforward correction

    % input limit active -> clamp feedforward, no feedback
    if U(k) + dk > Fmax || U(k) + dk < -Fmax
        dk = min(Fmax, max(-Fmax, U(k) + dk)) - U(k);
        Kk = zeros(1,4);
    end

    % cost-to-go coefficients at time k (uses S_{k+1}, s_{k+1})
    s = Q*X(:,k) + Ak'*s + Ak'*S*Bk*dk;
    S = Q + Ak'*S*(Ak - Bk*Kk);
    S = 0.5*(S + S');                   % keep symmetric (numerics)

    K(:,:,k) = Kk;
    d(k)     = dk;
end

end
