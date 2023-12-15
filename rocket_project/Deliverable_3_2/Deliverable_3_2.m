addpath(fullfile('..', 'src'));

close all
clear all
clc

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/20;

% Choose simulation time
Tf = 10;
rocket= Rocket(Ts); 
[xs,us] =rocket.trim(); %Compute steady−state for which 0=f(xs,us) 
sys=rocket.linearize(xs,us); %Linearizethenonlinearmodelabouttrimpoint

[sysx, sysy, sysz, sysroll]=rocket.decompose(sys,xs,us);

% Reference states
pos_ref = -4;
angle_ref = deg2rad(35);
%% sysx
Hx = 5; % Horizon length in seconds
mpc_x = MpcControl_x(sysx, Ts, Hx);
x0_x = [0;0;0;3];
u_x = mpc_x.get_u(x0_x, pos_ref);

[T, X_sub, U_sub] = rocket.simulate_f(sysx, x0_x, Tf, @mpc_x.get_u, pos_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sysx, xs, us, pos_ref);

%% sysy
Hy = 5; % Horizon length in seconds
mpc_y = MpcControl_y(sysy, Ts, Hy);
x0_y = [0;0;0;3];
u_y = mpc_y.get_u(x0_y, pos_ref);

[T, X_sub, U_sub] = rocket.simulate_f(sysy, x0_y, Tf, @mpc_y.get_u, pos_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sysy, xs, us, pos_ref);

%% Pdiff
Hroll = 5; % Horizon length in seconds
mpc_roll = MpcControl_roll(sysroll, Ts, Hroll);
x0_roll = [0;deg2rad(30)];
u_roll = mpc_roll.get_u(x0_roll, angle_ref);

[T, X_sub, U_sub] = rocket.simulate_f(sysroll, x0_roll, Tf, @mpc_roll.get_u, angle_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sysroll, xs, us, angle_ref);

%% Pavg
Hz = 5; % Horizon length in seconds
mpc_z = MpcControl_z(sysz, Ts, Hz);
x0_z = [0;3];
u_z = mpc_z.get_u(x0_z, pos_ref);

[T, X_sub, U_sub] = rocket.simulate_f(sysz, x0_z, Tf, @mpc_z.get_u, pos_ref);
ph = rocket.plotvis_sub(T, X_sub, U_sub, sysz, xs, us, pos_ref);