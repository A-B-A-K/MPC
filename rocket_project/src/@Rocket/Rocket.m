classdef Rocket
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rocket class
    %
    % Structure (Click VIEW > Collapse All for better overview)
    %
    % methods
    % - System
    % - Simulation
    % - Visualization
    %
    % methods (hidden)
    % - System
    % - Simulation
    % - Visualization
    %
    % methods (static)
    % - Helpers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        sys        % State, input, output names and units
        Ts         % Sampling time
        delay = 0; % Delay in samples for control signal to reach the actuators
        
        anim_rate = 0; % Animation rate for visualization, default = 0x.
        
        % Physical properties and constraints
        mass      = 1.7;
        mass_rate = 0;                                   % Change of mass per second and throttle^2
        g         = 9.81;                                % gravitational acceleration
        J         = diag([0.0644; 0.0644; 0.0128]);      % inertia tensor
        inv_J     = inv(diag([0.0644; 0.0644; 0.0128]));
        
        % Propeller model constants
        thrust_coeff = [0; 0.03; 0]; % experimentally identified
        torque_coeff = -0.1040;           % experimentally identified
        r_F          = [0; 0; -0.215];    % thruster position in body frame
    end
    properties (Constant)
        % System dimensions
        dimx = 12; % Length of state vector
        dimu = 4;  % Length of input vector
        
        % Input limits
        ubu = [deg2rad([ 15;  15]); 80;  20]; % upper bound for input: [rad, rad, %, %]
        lbu = [deg2rad([-15; -15]); 50; -20]; % lower bound for input: [rad, rad, %, %]
        
        % State limits for (linear) control
        ubx = [deg2rad([Inf Inf Inf, 10 10 Inf]), Inf Inf Inf, Inf Inf Inf]';
        lbx = -[deg2rad([Inf Inf Inf, 10 10 Inf]), Inf Inf Inf, Inf Inf Inf]';
        
        % Indices in the state vector
        indx = struct('omega', 1:3, 'phi', 4:6, 'vel', 7:9, 'pos', 10:12);
        indu = struct('d1', 1, 'd2', 2, 'Pavg', 3, 'Pdiff', 4);
    end
    properties (Hidden, Constant)
        % Plotting colors
        color = struct('meas', [0, 0.4470, 0.7410], 'ref', 'k');
    end
    
    methods
        %
        % Constructor
        %
        function obj = Rocket(Ts)
            
            obj.Ts = Ts;
            
            try % Test YALMIP installation
                sdpvar(2,1);
            catch
                error('Could not load YALMIP - check that it is installed properly and on the path.')
            end
            try % Test casadi installation
                import casadi.*
                casadi.SX.sym('x');
            catch
                error('Could not load casadi - check that it is installed properly and on the path.')
            end
            
            % Define system state, input, output names and inputs
            sys.StateName = { ...
                'wx', 'wy', 'wz', ...
                'alpha', 'beta', 'gamma',...
                'vx', 'vy', 'vz', ...
                'x','y','z'};
            sys.StateUnit = { ...
                'rad/s', 'rad/s', 'rad/s', ...
                'rad', 'rad', 'rad', ...
                'm/s', 'm/s', 'm/s', ...
                'm', 'm', 'm'};
            
            sys.InputName = {'d1', 'd2', 'Pavg', 'Pdiff'};
            sys.InputUnit = {'rad', 'rad', '%', '%'};
            
            sys.OutputName = sys.StateName;
            sys.OutputUnit = sys.StateUnit;
            
            obj.sys = sys;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % System
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Continuous rocket dynamics x_dot = f(x,u)
        % ----------------------------------------------------
        % x = [w; phi; v; p]
        % w     [rad/s] Angular velocity in body frame (b, {forward, left, up})
        % phi   [rad]   Euler angles expressing the attitude of body frame w.r.t. world frame (w, {East, North, Up})
        % v     [m/s]   Linear velocity in world frame
        % p     [m]     Position in world frame
        % ----------------------------------------------------
        % u = [d1; d2; P_avg; P_delta]
        % d1    [rad] Thrust vector deflection producing negativ angular velocity about body x axis
        % d2    [rad] Thrust vector deflection producing negativ angular velocity about body y axis
        % Pavg  [%]   Avg propeller throttle
        % Pdiff [%]   Difference propeller throttle producing negative angular velocity about body z axis
        %
        function [x_dot, y] = f(obj, x, u)
            
            if nargin < 3
                u = zeros(4, 1);
            end
            
            % Decompose state
            [w, phi, v, ~] = obj.parse_state(x);
            
            % Get transformation from body to world frame
            Twb = obj.eul2mat(phi);
            
            % Get force and moment in body frame from (differential) thrust
            [b_F, b_M] = obj.getForceAndMomentFromThrust(u);
            
            % State derivative components
            % Angular velocity in body frame
            w_dot = obj.inv_J * (b_M - cross(w, obj.J * w));
            
            % Attitude angles. Angular velocities in body frame -> rate of
            % change of Euler angles
            bet = phi(2); gam = phi(3);
            
            E_inv  = 1/cos(bet) * ...
                [cos(gam),         -sin(gam),          0;
                sin(gam)*cos(bet), cos(gam)*cos(bet),  0;
                -cos(gam)*sin(bet), sin(gam)*sin(bet), cos(bet)];
            
            phi_dot = E_inv * w;
            
            % Linear velocity and position in world frame
            v_dot = Twb * b_F/obj.mass - [0; 0; obj.g];
            p_dot = v;
            
            % Compose state derivative
            x_dot = [w_dot; phi_dot; v_dot; p_dot];
            
            % Output thrust force and moment for debugging
            y = [b_F; b_M];
        end
        
        %
        % Get force and moment vector from thrust in body reference frame
        %
        function [b_F, b_M] = getForceAndMomentFromThrust(obj, u)
            
            % Thrust and torque magnitude from avg and diff throttle command
            thrust = obj.g * obj.thrust_coeff' * [u(3)^2; u(3); 1];
            torque = obj.torque_coeff * obj.J(3,3) * u(4);
            
            % Direction of thrust and differential torque vector in body frame (unit vector)
            b_eF = [sin(u(2)); -sin(u(1))*cos(u(2)); cos(u(1))*cos(u(2))];
            
            % Resulting force and torque vector in body frame
            b_F = thrust * b_eF;
            b_M = torque * b_eF + cross(obj.r_F, b_F);
        end
        
        %
        % Compute trim point (steady flight, x_dot = 0)
        %
        function [xs, us] = trim(obj)
            
            % Since we know here that a trim point f(xs, us) = 0 exists,
            % we can simply solve an box-constrained optimization problem
            % with Matlab's fmincon solver:
            %                   min f(xs,us).^2
            %
            % Given that we have nonlinear dynamics, we will find the local
            % minimum next to the initial guess.
            %
            % Decision variable y = [x; u]
            
            MSE = @(v) v'*v;
            objective_fcn = @(y) MSE( obj.f(y(1:12), y(13:16)) );
            
            % Set bounds and initial guess
            ubx_ = [Inf Inf Inf, deg2rad([180 89 180]), Inf Inf Inf, Inf Inf Inf]';
            lbx_ = -ubx_;
            
            uby = [ubx_; obj.ubu(:)];
            lby = [lbx_; obj.lbu(:)];
            
            y_guess = zeros(16,1);
            
            % Setup solver
            opt = optimoptions('fmincon', 'Algorithm','sqp');
            opt.Display = 'off';
            
            % Solve
            [y, fval, exitflag] = fmincon(objective_fcn, y_guess, ...
                [], [], [], [], lby, uby, [], opt);
            
            if exitflag < 0 || fval > 1e-3
                error('Could not find trim condition');
            end
            xs = y(1:12);
            us = y(13:16);
            
            % Clean up numerical inaccuracies
            xs(abs(xs) < 5e-2) = 0;
            us(abs(us) < 1e-3) = 0;
        end
        
        %
        % Return the linearization of the system around the
        % equilibrium point (xs, us)
        %
        function linSys = linearize(obj, xs, us)
            if nargin < 3
                fprintf('No equilibrium given... trimming\n');
                [xs, us] = obj.trim();
            end
            
            % Create casadi symbolic variables
            x = casadi.SX.sym('x', 12);
            u = casadi.SX.sym('u', 4);
            f = obj.f(x, u);
            
            % Create symbolic casadi function for automatic differentiation
            % of A = df/dx, B = df/du. Evaluate and densify.
            A = casadi.Function('A', {x,u}, {jacobian(f, x)});
            A = full(A(xs, us));
            B = casadi.Function('B', {x,u}, {jacobian(f, u)});
            B = full(B(xs, us));
            
            % Clean up numerical inaccuracies
            A(abs(A) < 1e-5) = 0;
            B(abs(B) < 1e-5) = 0;
            
            % Create state space representation
            linSys = ss(A, B, eye(12), zeros(12,4));
            linSys.UserData.xs = xs;
            linSys.UserData.us = us;
            
            linSys.StateName = obj.sys.StateName;
            linSys.StateUnit = obj.sys.StateUnit;
            
            linSys.InputName = obj.sys.InputName;
            linSys.InputUnit = obj.sys.InputUnit;
            
            linSys.OutputName = obj.sys.OutputName;
            linSys.OutputUnit = obj.sys.OutputUnit;
        end
        
        %
        % Decompose the system into four systems around a hovering
        % equilibrium
        %
        function [sys_x, sys_y, sys_z, sys_roll] = decompose(obj, linSys, xs, us)
            
            % Split into four seperate systems
            I = obj.indx;
            
            % [wy beta vx x], [d2]
            idx = [I.omega(2) I.phi(2) I.vel(1) I.pos(1)];
            idu = 2;
            idy = I.pos(1);
            
            sys_x = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys x');
            
            % [wx alpha vy y], [d1]
            idx = [I.omega(1) I.phi(1) I.vel(2) I.pos(2)];
            idu = 1;
            idy = I.pos(2);
            
            sys_y = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys y');
            
            % [vz z], [Pavg]
            idx = [I.vel(3) I.pos(3)];
            idu = 3;
            idy = I.pos(3);
            
            sys_z = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys z');
            
            % [wz gamma], [Pdiff]
            idx = [I.omega(3) I.phi(3)];
            idu = 4;
            idy = I.phi(3);
            
            sys_roll = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys roll');
        end
        
        %
        % Decompose the system into four systems around a hovering
        % equilibrium
        %
        function mpc = merge_lin_controllers(obj, xs, us, mpc_x, mpc_y, mpc_z, mpc_r)
            
            % Get state indices
            linSys = obj.linearize(xs, us);
            [sys_x, sys_y, sys_z, sys_roll] = obj.decompose(linSys, xs, us);
            
            Iu = zeros(4, 1);
            
            idx_x = sys_x.UserData.idx;
            idu_x = sys_x.UserData.idu;
            Iu_x  = Iu; Iu_x(idu_x) = 1;
            
            idx_y = sys_y.UserData.idx;
            idu_y = sys_y.UserData.idu;
            Iu_y  = Iu; Iu_y(idu_y) = 1;
            
            idx_z = sys_z.UserData.idx;
            idu_z = sys_z.UserData.idu;
            Iu_z  = Iu; Iu_z(idu_z) = 1;
            
            idx_r = sys_roll.UserData.idx;
            idu_r = sys_roll.UserData.idu;
            Iu_r  = Iu; Iu_r(idu_r) = 1;
            
            % Define a local function so we can set a break point in the
            % function to evaluate sub-controllers independently
            function u = merged_get_u(z_, ref_)
                % If z_ is the state vector, 13:end = [], and mpc_z will be
                % evaluated with 0 disturbance
                u_x = mpc_x.get_u( z_(idx_x) - xs(idx_x),              ref_(1) );
                u_y = mpc_y.get_u( z_(idx_y) - xs(idx_y),              ref_(2) );
                u_z = mpc_z.get_u([z_(idx_z) - xs(idx_z); z_(13:end)], ref_(3) );
                u_r = mpc_r.get_u( z_(idx_r) - xs(idx_r),              ref_(4) );
                
                delta_u = Iu_x .* u_x + Iu_y .* u_y + Iu_z .* u_z + Iu_r .* u_r;
                u = us + delta_u;
            end
            
            % Define a local function so we can set a break point in the
            % function to evaluate sub-controllers independently
            function [u, T_opt, X_opt, U_opt] = merged_TXU_opt(z_, ref_)
                % If z_ is the state vector, 13:end = [], and mpc_z will be
                % evaluated with 0 disturbance
                [u_x, T_opt_x, X_opt_x, U_opt_x] = mpc_x.get_u( z_(idx_x) - xs(idx_x),              ref_(1) );
                [u_y, T_opt_y, X_opt_y, U_opt_y] = mpc_y.get_u( z_(idx_y) - xs(idx_y),              ref_(2) );
                [u_z, T_opt_z, X_opt_z, U_opt_z] = mpc_z.get_u([z_(idx_z) - xs(idx_z); z_(13:end)], ref_(3) );
                [u_r, T_opt_r, X_opt_r, U_opt_r] = mpc_r.get_u( z_(idx_r) - xs(idx_r),              ref_(4) );
                
                delta_u = Iu_x .* u_x + Iu_y .* u_y + Iu_z .* u_z + Iu_r .* u_r;
                u = us + delta_u;
                
                % Fuse different horizon-length trajectories
                Nx_xyzr = [size(X_opt_x, 2); size(X_opt_y, 2); size(X_opt_z, 2); size(X_opt_r, 2)];
                delta_X_opt = nan(obj.dimx, max(Nx_xyzr));
                delta_X_opt(idx_x,1:Nx_xyzr(1)) = X_opt_x;
                delta_X_opt(idx_y,1:Nx_xyzr(2)) = X_opt_y;
                delta_X_opt(idx_z,1:Nx_xyzr(3)) = X_opt_z;
                delta_X_opt(idx_r,1:Nx_xyzr(4)) = X_opt_r;
                X_opt = xs + delta_X_opt;
                                
                Nu_xyzr = [size(U_opt_x, 2); size(U_opt_y, 2); size(U_opt_z, 2); size(U_opt_r, 2)];
                delta_U_opt = nan(obj.dimu, max(Nu_xyzr));
                delta_U_opt(idu_x,1:Nu_xyzr(1)) = U_opt_x;
                delta_U_opt(idu_y,1:Nu_xyzr(2)) = U_opt_y;
                delta_U_opt(idu_z,1:Nu_xyzr(3)) = U_opt_z;
                delta_U_opt(idu_r,1:Nu_xyzr(4)) = U_opt_r;
                U_opt = us + delta_U_opt;
                
                % Pick longest time grid
                T_opt = T_opt_x;
                if T_opt_y(end) > T_opt(end), T_opt = T_opt_y; end
                if T_opt_z(end) > T_opt(end), T_opt = T_opt_z; end
                if T_opt_r(end) > T_opt(end), T_opt = T_opt_r; end
                
                if (max(Nx_xyzr) ~= min(Nx_xyzr) || ...
                    max(Nu_xyzr) ~= min(Nu_xyzr))
                    warning('Merged trajectories of different lengths. Some states might not be visualized correctly.');
                end
            end
            
            function [u, T_opt, X_opt, U_opt] = get_u_wrapper(z_, ref_)
                % If only u requested, use lightweight function to avoid
                % unnecessary evaluations of T_opt, X_opt, U_opt
                if nargout >= 2
                    [u, T_opt, X_opt, U_opt] = merged_TXU_opt(z_, ref_);
                else
                    u = merged_get_u(z_, ref_);
                end
            end
            
            % Return function handle
            mpc.get_u = @get_u_wrapper;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Simulate any system f in steps of Ts using an open-loop input
        % trajectory U, a single input u, or a control law @(x_, ref_)
        %
        [T, X, U, Ref, Z_hat] = simulate_f(obj, f, x0, Tf, U, ref, esti)
        
        %
        % Simulate the nonlinear model in steps of obj.Ts using an open-loop
        % input trajectory U, a single input u, or a control law @(x_, ref_)
        %
        [T, X, U, Ref, Z_hat] = simulate(obj, x0, Tf, U, ref, esti)
        
        %
        % Simulate the nonlinear model in steps of obj.Ts using an open-loop
        % input trajectory U, a single input u, or a control law @(x_, ref_).
        % Use hardcoded estimator in z direction.
        %
        [T, X, U, Ref, Z_hat] = simulate_est_z(obj, x0, Tf, U, ref, mpc_z, sys_z)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Visualization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Visualize the trajectory (without plots)
        %
        ph = vis(obj, T, X, U, ax)
        
        %
        % Plot & visualize full trajectory
        %
        ph = plotvis(obj, T, X, U, Ref)
        
        %
        % Plot and visualize sub-system trajectory (fill up other states
        % with (xs, us) for 3D visualization
        %
        ph = plotvis_sub(obj, T, subX, subU, subSys, xs, us, ref)
        
    end
    
    methods (Hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % System
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Continuous rocket dynamics x_dot = f(x, u) along 1D z direction
        %
        function x_dot = f_z(obj, x, u)
            
            if nargin < 3
                u = 1;
            end
            
            % Decompose state
            vz = x(1);
            
            thrust = obj.g * obj.thrust_coeff' * [u^2; u; 1];
            
            % State derivative components
            % Vertical speed and height
            vz_dot = thrust/obj.mass - obj.g;
            z_dot = vz;
            
            % Compose state derivative
            x_dot = [vz_dot; z_dot];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulation (hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Simulate the system xdot = f(x,u) from x0 forward for dt seconds
        %
        function xp = simulate_step(~, f, x0, u, dt)
            % Integrate forward to next time step
            [~, xout] = ode45( @(t_, x_) f(x_, u), [0, dt], x0 );
            xp = xout(end,:)';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Visualization (hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a graphical representation of the rocket
        %
        [h_rocket, h_thrust] = create_rocket_transformobj(obj, ax)
        
        %
        % Plot center of gravity point(s) into ax
        %
        visualize_pos(obj, pos, ax)
        
        %
        % Visualize the rocket at a given state and input
        %
        ph = visualize_point(obj, x, u, ax, h_rocket, h_thrust)
        
        %
        % Plot content into predefined axes
        %
        plot_into_axes(obj, ax, ux_id, sys, T, X, U, bx, bu, X_ref)
        
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helpers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Split the state into its parts, in: State trajectory (nx x time)
        %
        function [omega, phi, vel, pos] = parse_state(X)
            if nargout >= 1, omega = X(Rocket.indx.omega, :); end
            if nargout >= 2, phi   = X(Rocket.indx.phi,   :); end
            if nargout >= 3, vel   = X(Rocket.indx.vel,   :); end
            if nargout >= 4, pos   = X(Rocket.indx.pos,   :); end
        end
        
        %
        % Obtain transformation matrix from attitude angles
        %
        function Twb = eul2mat(eul)
            alp = eul(1); % alpha
            bet = eul(2); % beta
            gam = eul(3); % gamma
            
            % T1 is elementary rotation about x axis
            function T = T1(a)
                T = [1 0 0;
                    0 cos(a) sin(a);
                    0 -sin(a) cos(a)];
            end
            % T2 is elementary rotation about y axis
            function T = T2(a)
                T = [cos(a) 0 -sin(a);
                    0 1 0;
                    sin(a) 0 cos(a)];
            end
            % T3 is elementary rotation about z axis
            function T = T3(a)
                T = [cos(a) sin(a) 0;
                    -sin(a) cos(a) 0;
                    0 0 1];
            end
            % %             % Detailed derivation
            % %
            % %             % Going from world to body frame:
            % %             % 1. Rotate by alpha about world x axis:
            % %             T_bw1 = T1(alp);
            % %             % 2. Rotate by beta about resulting y axis;
            % %             T_bw2 = T2(bet) * T_bw1;
            % %             % 3. Rotate by gamma about resluting z axis;
            % %             T_bw = T3(gam) * T_bw2;
            % %
            % %             % In one step:
            % %             T_bw = T3(gam) * T2(bet) * T1(alp);
            % %
            % %             % T_bw transforms a w_vect in body frame (b_vect):
            % %             % b_vect = T_bw * w_vect
            % %             % We need the opposite, so we transpose the matrix:
            % %             T_wb = T_bw';
            % %
            % %             % Alternatively, we can directly construct this matrix by going
            % %             % backwards from body to world frame (same logic, inverse
            % %             % angles):
            Twb = T1(-alp) * T2(-bet) * T3(-gam);
        end
        
        %
        % Get safe outer figure size for specific screen
        %
        function position = get_position_on_screen(targetSize, max_screen_usage)
            
            if nargin < 2
                max_screen_usage = 0.9;
            end
            
            figRatio = targetSize(1) / targetSize(2);
            
            screenSize = get(0, 'ScreenSize'); screenSize = screenSize(3:4);
            availSize = max_screen_usage * screenSize;
            horClip = [availSize(1), availSize(1) / figRatio];
            verClip = [figRatio * availSize(2), availSize(2)];
            minClip = min(horClip, verClip);
            figSize = round( min(minClip, targetSize) );
            
            % Determine figure offset (window/menu bar)
            f = figure('Visible', 'off');
            sizeOffset = f.OuterPosition(3:4) - f.Position(3:4);
            close(f);
            
            outerPos = screenSize - figSize - sizeOffset - [0 50]; % 50 height safety for Windows
            position = [outerPos, figSize];
        end
        
        %
        % Create pose axes (3D, equal, hold, grid, labels)
        %
        ph = create_pose_axes(figNo_or_ax)
        
        %
        % Create sub-system from indices
        %
        function sub_sys = parse_system(sys, xs, us, idx, idu, idy, name)
            
            [A, B, C, ~] = ssdata(sys);
            
            sub_sys = ss(A(idx,idx), B(idx,idu), C(idy,idx), 0);
            
            sub_sys.UserData.idx = idx;
            sub_sys.UserData.idu = idu;
            sub_sys.UserData.idy = idy;
            
            sub_sys.UserData.xs = xs(idx);
            sub_sys.UserData.us = us(idu);
            
            sub_sys.Name = name;
            sub_sys.StateName = sys.StateName(idx);
            sub_sys.StateUnit = sys.StateUnit(idx);
            
            sub_sys.InputName = sys.InputName(idu);
            sub_sys.InputUnit = sys.InputUnit(idu);
            
            sub_sys.OutputName = sys.OutputName(idy);
            sub_sys.OutputUnit = sys.OutputUnit(idy);
        end
       
    end
    
end