addpath(fullfile('..', 'src'));

close all
clear all
clc

%% TODO: This file should produce all the plots for the deliverable

% Define subsystems
Ts = 1/20;

% Choose simulation time
Tf = 10;
rocket= Rocket(Ts); 
[xs,us] =rocket.trim(); %Compute steady−state for which 0=f(xs,us) 
sys=rocket.linearize(xs,us); %Linearize the non linear model about trim point

[sysx, sysy, sysz, sysroll]=rocket.decompose(sys,xs,us);

H = 5;
mpc_x = MpcControl_x(sysx, Ts, H);
mpc_y = MpcControl_y(sysy, Ts, H);
mpc_z = MpcControl_z(sysz, Ts, H);
mpc_roll = MpcControl_roll(sysroll, Ts, H);

% Merge four sub−system controllers into one full−system controller
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
x0 = [zeros(1, 9), 1 0 3]';
ref = [1.2, 0, 3, 0]';

% Manipulate mass for simulation
rocket.mass = 2.13;

% Using state estimates of observer
% [T, X, U, Ref, Z_hat] = rocket.simulate_est_z(x0, Tf, @mpc.get_u, ref, mpc_z, sysz);
% rocket.anim_rate = 5;
% ph = rocket.plotvis(T, X, U, Ref);
% % 
% % Not using state estimates from observer
rocket.mass = 2.13;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @mpc.get_u, ref);
rocket.anim_rate = 5;
ph = rocket.plotvis(T, X, U, Ref);