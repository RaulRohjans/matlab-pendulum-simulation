%% Physical parameters 
%%%%% These are the parameters used by the Simulink model and they
%%%%% correspond to the properties of the real physical system

p.m = 1.0;        % kg
p.l = 0.5;        % m, distance pivot -> point mass
p.b = 0.05;       % N*m/(rad/s)
g = 9.81;       % m/s^2

p.theta0 = 35;  %degrees
p.omega0 = 0;   %rad/s

%% Model parameters 
%%%%%%% These are the parameters you can use for your model. Currently,
%%%%%%% they are set equal to the values of the physical system but you can
%%%%%%% change and observe how your controller and observer behave under model
%%%%%%% uncertainty

m.m = p.m;        % kg
m.l = p.l;        % m, distance pivot -> point mass
m.b = p.b;       % N*m/(rad/s)
m.J = m.m*m.l^2;


%% Sampling
%%%%%% Controller works with measurements obtained at regular intervals.
%%%%%% The control inputs are in discrete time as well. They are then
%%%%%% applied to the system in zero-order-hold manner. 
%%%%%% Sampling period: every Ts seconds
Ts = 0.02; 

%% Actuator limits 
umin = ...
umax = ...

%% Linear model
%%%%%% Here you should set the matrices Ac and Bc which correspond to
%%%%%% the linearized system 

Ac = ...
Bc = ...

%% Discrete-time model
%%%%% Disretize the linear apprroximation; you can use MATLAB functions
%%%%% from control systtem toolbox such as c2d or compute matrix
%%%%% exponentials

Ad=...
Bd=...

%% LQR
%%%%% Design LQR 

K=...


%% Observer
%%%% Design an observer using pole placement
C = [1 0];
xhat0 = ... 

L = ...

%% Simulate the actuated pendulum
%%% The following command runs the simulation.
%%% You should see an animation of the pendulum and out is a struct that
%%% contains the outputs from the simulation (for example, measured and
%%% estimated values, control input, or any other thing you add)

out=sim("actuated_pendulum")

%%%%%% 
%%% Add relevant plots etc. 


