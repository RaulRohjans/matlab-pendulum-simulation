function [A,B] = linearize_trajectory(X,U,model,Ts)
%LINEARIZE_TRAJECTORY Linearize the discrete nonlinear dynamics numerically.
% A(:,:,k) and B(:,:,k) are central-difference Jacobians of discrete_step
% about the accepted nominal pair X(:,k), U(k).

N = numel(U);
assert(isequal(size(X), [4, N+1]), ...
    'X must have size 4-by-(numel(U)+1).');
A = zeros(4,4,N);
B = zeros(4,1,N);

for k = 1:N
    x = X(:,k);
    u = U(k);

    for j = 1:4
        h = 1e-5*max(1,abs(x(j)));
        e = zeros(4,1);
        e(j) = h;

        A(:,j,k) = (discrete_step(x+e,u,model,Ts) ...
            - discrete_step(x-e,u,model,Ts))/(2*h);
    end

    hu = 1e-5*max(1,abs(u));
    B(:,:,k) = (discrete_step(x,u+hu,model,Ts) ...
               - discrete_step(x,u-hu,model,Ts))/(2*hu);
end
end
