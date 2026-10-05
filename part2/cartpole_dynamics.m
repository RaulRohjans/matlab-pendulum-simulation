function dx = cartpole_dynamics(x,u,p)
%CARTPOLE_DYNAMICS Continuous nonlinear cart-pole dynamics.
% State: x = [cart position; cart velocity; pole angle; angular velocity].
% The pole angle is zero when upright and positive counterclockwise.

x = x(:);
u = u(1);
v = x(2);
theta = x(3);
omega = x(4);

H = [p.M+p.m,  -p.m*p.ell*cos(theta);
     -p.m*p.ell*cos(theta),   p.m*p.ell^2];

rhs = [u-p.bc*v-p.m*p.ell*sin(theta)*omega^2;
       p.m*p.g*p.ell*sin(theta)-p.bp*omega];

acc = H\rhs;
dx = [v; acc(1); omega; acc(2)];
end
