function xn = discrete_step(x,u,p,Ts)
%DISCRETE_STEP Propagate the nonlinear dynamics for one sample with RK4.
% The applied force is constant throughout the sample (zero-order hold).

x = x(:);
a = cartpole_dynamics(x, u, p);
b = cartpole_dynamics(x + Ts*a/2, u, p);
c = cartpole_dynamics(x + Ts*b/2, u, p);
d = cartpole_dynamics(x + Ts*c, u, p);

xn = x + Ts*(a+2*b+2*c+d)/6;
end
