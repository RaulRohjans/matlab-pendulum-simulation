function [A,B] = linearize_trajectory(X,U,model,Ts)
% Numerical derivatives for the discrete transition

N = numel(U);
A = zeros(4,4,N);
B = zeros(4,1,N);
h = 1e-5;

for k = 1:N
    x = X(:,k);
    u = U(k);

    for j = 1:4
        e = zeros(4,1);
        e(j) = h;

        A(:,j,k) = (discrete_step(x+e,u,model,Ts) ...
-discrete_step(x-e,u,model,Ts))/(2*h);
    end

    B(:,:,k) = (discrete_step(x,u+h,model,Ts) ...
               -discrete_step(x,u-h,model,Ts))/(2*h);
end
end