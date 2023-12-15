addpath(fullfile('..', 'src'));

%close all
%clear all
%clc

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/20;

% Choose simulation time
Tf = 8;
rocket= Rocket(Ts); 
[xs,us] =rocket.trim(); %Compute steady−state for which 0=f(xs,us) 
sys=rocket.linearize(xs,us); %Linearizethenonlinearmodelabouttrimpoint

[sysx, sysy, sysz, sysroll]=rocket.decompose(sys,xs,us);

% Design MPC controller
Hx = 5; % Horizon length in seconds
mpc_x = MpcControl_x(sysx, Ts, Hx);

% Get control input (x is the index of the subsystem here)
x0_x = [0;0;0;3];
u_x = mpc_x.get_u(x0_x);

% Set initial state stationary/not rotating at position x=3.
[T, X_sub, U_sub] = rocket.simulate_f(sysx, x0_x, Tf, @mpc_x.get_u, 0);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sysx, xs, us);
