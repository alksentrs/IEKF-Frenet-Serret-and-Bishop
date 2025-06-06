classdef bishop_ekf < handle
    properties
        % State variables
        x_hat     % Position estimate (3x1)
        R_hat     % Rotation matrix estimate (3x3)
        kappa1_hat % First Bishop curvature estimate (scalar)
        kappa2_hat % Second Bishop curvature estimate (scalar)
        u_hat     % Tangential velocity estimate (scalar)
        a_hat     % Tangential acceleration estimate (scalar)
        
        % Covariance matrix
        P         % Error covariance matrix (10x10)
        
        % Process noise parameters
        Q_x       % Position process noise (3x3)
        Q_omega   % Angular velocity process noise (3x3)
        Q_kappa1  % First curvature process noise (scalar)
        Q_kappa2  % Second curvature process noise (scalar)
        Q_u       % Tangential velocity process noise (scalar)
        Q_a       % Tangential acceleration process noise (scalar)
        
        % Measurement noise
        R         % Measurement noise covariance (3x3)
        
        % Time step
        dt        % Time step for integration
    end
    
    methods
        function obj = bishop_ekf(dt, Q_x, Q_omega, Q_kappa1, Q_kappa2, Q_u, Q_a, R)
            % Constructor
            obj.dt = dt;
            obj.Q_x = Q_x;
            obj.Q_omega = Q_omega;
            obj.Q_kappa1 = Q_kappa1;
            obj.Q_kappa2 = Q_kappa2;
            obj.Q_u = Q_u;
            obj.Q_a = Q_a;
            obj.R = R;
            
            % Initialize state
            obj.x_hat = zeros(3,1);
            obj.R_hat = eye(3);
            obj.kappa1_hat = 0.08;  % Initial first curvature
            obj.kappa2_hat = 0.03;  % Initial second curvature
            obj.u_hat = 8.0;        % Initial tangential velocity
            obj.a_hat = 0.01;        % Initial tangential acceleration
            
            % Initialize covariance (10x10 for all state variables)
            obj.P = eye(10);
        end
        
        function [x_hat, R_hat] = predict(obj)
            % Prediction step
            % Update position
            v_body = [obj.u_hat; 0; 0];
            obj.x_hat = obj.x_hat + obj.R_hat * v_body * obj.dt;
            
            % Update rotation matrix using Bishop frame formulas
            omega = [0; -obj.kappa2_hat; obj.kappa1_hat];  % Bishop Darboux vector
            Omega = skew(omega);
            obj.R_hat = obj.R_hat * expm(Omega * obj.dt);
            
            % Update velocity and acceleration
            obj.u_hat = obj.u_hat + obj.a_hat * obj.dt;
            
            % Update covariance
            F = obj.compute_jacobian();
            Q = blkdiag(obj.Q_omega, obj.Q_x, obj.Q_kappa1, obj.Q_kappa2, obj.Q_u, obj.Q_a);
            obj.P = F * obj.P * F' + Q * obj.dt;
            
            x_hat = obj.x_hat;
            R_hat = obj.R_hat;
        end
        
        function update(obj, y)
            % Update step
            % Compute innovation
            innovation = y - obj.x_hat;
            
            % Compute measurement Jacobian (3x10)
            H = [zeros(3,3), eye(3), zeros(3,4)];
            
            % Compute Kalman gain
            S = H * obj.P * H' + obj.R;
            K = obj.P * H' / S;
            
            % Update state
            delta = K * innovation;
            delta_omega = delta(1:3);
            delta_x = delta(4:6);
            delta_kappa1 = delta(7);
            delta_kappa2 = delta(8);
            delta_u = delta(9);
            delta_a = delta(10);
            
            % Update rotation matrix
            obj.R_hat = obj.R_hat * expm(skew(delta_omega));
            
            % Update position using B(delta_omega)
            Bmat = obj.B(delta_omega);
            obj.x_hat = obj.x_hat + obj.R_hat * Bmat * delta_x;
            
            % Update curvatures, velocity and acceleration
            obj.kappa1_hat = obj.kappa1_hat + delta_kappa1;
            obj.kappa2_hat = obj.kappa2_hat + delta_kappa2;
            obj.u_hat = obj.u_hat + delta_u;
            obj.a_hat = obj.a_hat + delta_a;
            
            % Update covariance
            obj.P = (eye(10) - K * H) * obj.P;
        end
        
        function F = compute_jacobian(obj)
            % Compute system Jacobian matrix (10x10) for Bishop frame
            F = zeros(10,10);
            k1 = obj.kappa1_hat;
            k2 = obj.kappa2_hat;
            u = obj.u_hat;
            
            % The order is: [eta^R(3), eta^x(3), eta^kappa1, eta^kappa2, eta^u, eta^a]
            % Row 1 (first row of rotation error)
            F(1,2) = -k1;
            F(1,3) = -k2;
            
            % Row 2 (second row of rotation error)
            F(2,1) = k1;
            F(2,8) = 1;
            
            % Row 3 (third row of rotation error)
            F(3,1) = k2;
            F(3,7) = -1;
            
            % Row 4 (first row of position error)
            F(4,5) = -k1;
            F(4,6) = -k2;
            F(4,9) = -1;
            
            % Row 5 (second row of position error)
            F(5,3) = -u;
            F(5,4) = k1;
            
            % Row 6 (third row of position error)
            F(6,2) = u;
            F(6,4) = k2;
            
            % Row 9 (Velocity)
            F(9,10) = -1;
            
            % Row 10 (Acceleration)
            % F(10,:) = 0 by default
        end
        
        function Bmat = B(~, delta_omega)
            % Compute the B matrix for the update step
            theta = norm(delta_omega);
            if theta < 1e-8
                Bmat = eye(3);
            else
                wx = skew(delta_omega);
                Bmat = eye(3) + (1 - cos(theta)) / theta^2 * wx + (theta - sin(theta)) / theta^3 * (wx^2);
            end
        end
    end
end 