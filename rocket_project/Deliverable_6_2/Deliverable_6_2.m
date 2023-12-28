addpath(fullfile('..', 'src'));

%close all
%clear all
%clc


Ts = 1/40;
rocket = Rocket(Ts);

H = 2;     % Horizon length in seconds
delay = 0;
nmpc = NmpcControl (rocket, H, delay);
x0 = zeros(12,1);
ref = [0.5,0,1,deg2rad(65)]';


% Close-loop simuation
Tf = 2,5;
rocket.mass = 1.75;
rocket.delay = 4; % 0 if not specified
rocket.anim_rate = 1;
[T,X,U,Ref] = rocket.simulate(x0,Tf, @nmpc.get_u,ref);
ph = rocket.plotvis(T, X, U, Ref);

