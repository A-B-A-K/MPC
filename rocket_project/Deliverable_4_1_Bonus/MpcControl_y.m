classdef MpcControl_y < MpcControlBase
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N_segs = ceil(H/Ts); % Horizon steps
            N = N_segs + 1;      % Last index in 1-based Matlab indexing

            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            % System Dynamics
            A = mpc.A;
            B = mpc.B;

            % Cost matrices
            Q = [1.5 0 0 0;
                0 0.8 0 0;
                0 0 0.1 0;
                0 0 0 0.1];

            R = 1;
            Q_slack = 1 * eye(4); 


            % Constraints
            % u in U = { u | Mu <= m }
            M = [1;-1]; 
            m = [deg2rad(15); deg2rad(15)];
            % x in X = { x | Fx <= f }
            F = [0 1 0 0; 0 0 0 0;0 -1 0 0; 0 0 0 0]; 
            f = [deg2rad(10); 0; deg2rad(10); 0];
            slack = sdpvar(4, N-1);

            % Compute LQR controller for unconstrained system
            [K,Qf,~] = dlqr(A,B,Q,R);

            % MATLAB defines K as -K, so invert its signal
            K = -K; 


            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            obj = 0;
            con = [];
            for i = 1:N-1
                con = [con, X(:,i+1) == A*X(:,i) + B*U(:,i)];
                % Soften state constraints with slack
                con = [con, F*(X(:,i)-x_ref) <= f + slack(:,i)];
                con = [con, M*(U(:,i)-u_ref) <= m];
                obj = obj + (X(:,i)-x_ref)'*Q*(X(:,i)-x_ref) + (U(:,i)-u_ref)'*R*(U(:,i)-u_ref);
                % Penalize slack variable use
                obj = obj + slack(:,i)' * Q_slack * slack(:,i);
            end
            % con = con + (Ff*X(:,N) <= ff);
            obj = obj + (X(:,N)-x_ref)'*Qf*(X(:,N)-x_ref);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref}, {U(:,1), X, U});
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            A = mpc.A;
            B = mpc.B;
            C = mpc.C;

            % Rewrite constraints
            % u in U = { u | Mu <= m }
            M = [1;-1]; 
            m = [deg2rad(15); deg2rad(15)];

            % x in X = { x | Fx <= f }
            F = [0 1 0 0 ; 0 -1 0 0]; 
            f = [deg2rad(10); deg2rad(10)];

            obj = us^2;
            con = [xs == A*xs+B*us, ref == C*xs, M*us <= m, F*xs <= f];
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
