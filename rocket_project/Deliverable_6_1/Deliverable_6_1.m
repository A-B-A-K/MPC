addpath(fullfile('..', 'src'));

%close all
%clear all
%clc


Ts = 1/20;
rocket = Rocket(Ts);

H = 5;     % Horizon length in seconds
nmpc = NmpcControl (rocket, H);

% MPC reference with default maximum roll = 15 deg
ref = @(t_,x_) ref_TVC(t_);

% MPC reference with default maximum roll = 50 deg
%roll_max = deg2rad(50);
%ref = @(t_,x_) ref_TVC(t_, roll_max);

% Evaluate once and plot optimal open-loop trajectory,
% pad last input to get consistent size with time and state
x0 = [0;0;0;0;0;0;0;0;0;0;0;0];
[u, T_opt, X_opt, U_opt] = nmpc.get_u(x0, ref(1)');
U_opt(:,end+1) = nan;
ph = rocket.plotvis(T_opt, X_opt, U_opt, ref(1)');

% Close-loop simuation
Tf = 15;
rocket.anim_rate = 1;
[T, X, U, Ref] = rocket.simulate(x0, Tf, @nmpc.get_u, ref);
ph = rocket.plotvis(T, X, U, Ref);

