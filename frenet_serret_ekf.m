classdef frenet_serret_ekf < handle
    properties
        % State variables
        x_hat    % Position estimate (3x1)
        R_hat    % Rotation matrix estimate (3x3)
        gamma_hat % Curvature estimate (scalar)
        tau_hat   % Torsion estimate (scalar)
        u_hat     % Tangential velocity estimate (scalar)
        
        % Covariance matrix
        P        % Error covariance matrix (9x9)
        
        % Process noise parameters
        Q_x      % Position process noise (3x3)
        Q_omega  % Angular velocity process noise (3x3)
        Q_gamma  % Curvature process noise (scalar)
        Q_tau    % Torsion process noise (scalar)
        Q_u      % Tangential velocity process noise (scalar)
        
        % Measurement noise
        R        % Measurement noise covariance (3x3)
        
        % Time step
        dt       % Time step for integration
    end
    
    methods
        function obj = frenet_serret_ekf(dt, Q_x, Q_omega, Q_gamma, Q_tau, Q_u, R)
            % Constructor
            obj.dt = dt;
            obj.Q_x = Q_x;
            obj.Q_omega = Q_omega;
            obj.Q_gamma = Q_gamma;
            obj.Q_tau = Q_tau;
            obj.Q_u = Q_u;
            obj.R = R;
            
            % Initialize state
            obj.x_hat = zeros(3,1);
            obj.R_hat = eye(3);
            obj.gamma_hat = 0.08;
            obj.tau_hat = 0.03;
            obj.u_hat = 8.;
            
            % Initialize covariance (9x9 for all state variables)
            obj.P = eye(9);
        end
        
        function [x_hat, R_hat] = predict(obj)
            % Prediction step
            % Update position
            v_body = [obj.u_hat; 0; 0];
            obj.x_hat = obj.x_hat + obj.R_hat * v_body * obj.dt;
            
            % Update rotation matrix using Frenet-Serret formulas
            omega = [obj.tau_hat; 0; obj.gamma_hat];
            Omega = skew(omega);
            obj.R_hat = obj.R_hat * expm(Omega * obj.dt);
            
            % Update covariance
            F = obj.compute_jacobian();
            Q = blkdiag(obj.Q_omega, obj.Q_x, obj.Q_gamma, obj.Q_tau, obj.Q_u);
            obj.P = F * obj.P * F' + Q * obj.dt;
            
            x_hat = obj.x_hat;
            R_hat = obj.R_hat;
        end
        
        function update(obj, y)
            % Update step
            % Compute innovation
            innovation = y - obj.x_hat;
            
            % Compute measurement Jacobian (3x9)
            H = [zeros(3,3), eye(3), zeros(3,3)];
            
            % Compute Kalman gain
            S = H * obj.P * H' + obj.R;
            K = obj.P * H' / S;
            
            % Update state
            delta = K * innovation;
            delta_omega = delta(1:3);
            delta_x = delta(4:6);
            delta_gamma = delta(7);
            delta_tau = delta(8);
            delta_u = delta(9);
            
            % Update rotation matrix
            obj.R_hat = obj.R_hat * expm(skew(delta_omega));
            
            % Update position using B(delta_omega)
            Bmat = obj.B(delta_omega);
            obj.x_hat = obj.x_hat + obj.R_hat * Bmat * delta_x;
            
            % Update curvature, torsion, and tangential velocity
            obj.gamma_hat = obj.gamma_hat + delta_gamma;
            obj.tau_hat = obj.tau_hat + delta_tau;
            obj.u_hat = obj.u_hat + delta_u;
            
            % Update covariance
            obj.P = (eye(9) - K * H) * obj.P;
        end
        
        function F = compute_jacobian(obj)
            % Compute system Jacobian matrix (9x9) as in the paper
            F = zeros(9,9);
            g = obj.gamma_hat;
            t = obj.tau_hat;
            u = obj.u_hat;
            % The order is: [eta^R(3), eta^x(3), eta^gamma, eta^tau, eta^u]
            % Fill in the matrix as in the paper
            % Row 1
            F(1,2) = -g;
            F(1,8) = -1;
            % Row 2
            F(2,1) = g;
            F(2,3) = -t;
            % Row 3
            F(3,2) = t;
            F(3,7) = -1;
            % Row 4
            F(4,5) = -g;
            F(4,9) = -1;
            % Row 5
            F(5,3) = -u;
            F(5,4) = g;
            F(5,6) = -t;
            % Row 6
            F(6,2) = u;
            F(6,5) = t;
            % The rest are zeros (for gamma, tau, u)
            % F(7:9, :) = 0 by default
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