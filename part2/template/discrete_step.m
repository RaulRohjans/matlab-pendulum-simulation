function xn = discrete_step(x,u,p,Ts)
% Propagate system dynamics for Ts using 4th order Runge-Kutta

a = cartpole_dynamics(x,u,p);
b = cartpole_dynamics(x+Ts*a/2,u,p);
c = cartpole_dynamics(x+Ts*b/2,u,p);
d = cartpole_dynamics(x+Ts*c,u,p);

xn = x + Ts*(a+2*b+2*c+d)/6;
end