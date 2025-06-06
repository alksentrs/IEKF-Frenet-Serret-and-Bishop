% Test script for Frenet-Serret EKF implementation

% Parameters
dt = 0.05;  % Reduced time step for better observability
T = 10;     % Simulation time
N = T/dt;   % Number of steps

% Process noise parameters
Q_x = diag([0.01; 0.01; .02]);     % Position process noise
Q_omega = 0.1 * eye(3);  % Angular velocity process noise
Q_gamma = .1;          % Increased curvature process noise
Q_tau = 1.;            % Increased torsion process noise
Q_u = .1;              % Increased tangential velocity process noise (scalar)

% Measurement noise
R = diag([0.1; 0.1; .2]);       % Reduced measurement noise for better observability

% Initialize EKF
ekf = frenet_serret_ekf(dt, Q_x, Q_omega, Q_gamma, Q_tau, Q_u, R);

% Initialize true state
x_true = zeros(3,1);
R_true = eye(3);
gamma_true = 0.1;  % Constant curvature
tau_true = 0.05;   % Constant torsion
u_true = 10;        % Tangential velocity (scalar)

% Storage for plotting
x_est = zeros(3,N);
x_true_storage = zeros(3,N);
y_storage = zeros(3,N);  % Store measurements
gamma_est = zeros(1,N);  % Store curvature estimates
tau_est = zeros(1,N);    % Store torsion estimates
u_est = zeros(1,N);      % Store tangential velocity estimates
t = 0:dt:T-dt;

% Main loop
for k = 1:N
    % True state propagation
    omega = [tau_true; 0; gamma_true];
    Omega = skew(omega);
    R_true = R_true * expm(Omega * dt);
    x_true = x_true + R_true * [u_true; 0; 0] * dt;
    
    % Generate measurement
    y = x_true + sqrt(R) * randn(3,1);
    
    % EKF prediction
    [x_pred, R_pred] = ekf.predict();
    
    % EKF update
    ekf.update(y);
    
    % Store results
    x_est(:,k) = ekf.x_hat;
    x_true_storage(:,k) = x_true;
    y_storage(:,k) = y;
    gamma_est(k) = ekf.gamma_hat;
    tau_est(k) = ekf.tau_hat;
    u_est(k) = ekf.u_hat;
end

% Plot position results
figure('Name', 'Position Estimation');
subplot(3,1,1);
plot(t, x_est(1,:), 'b-', t, x_true_storage(1,:), 'r--', t, y_storage(1,:), 'k.');
ylabel('x');
legend('Estimated', 'True', 'Measurements');
grid on;

subplot(3,1,2);
plot(t, x_est(2,:), 'b-', t, x_true_storage(2,:), 'r--', t, y_storage(2,:), 'k.');
ylabel('y');
grid on;

subplot(3,1,3);
plot(t, x_est(3,:), 'b-', t, x_true_storage(3,:), 'r--', t, y_storage(3,:), 'k.');
ylabel('z');
xlabel('Time (s)');
grid on;

% 3D plot
figure('Name', '3D Trajectory');
plot3(x_est(1,:), x_est(2,:), x_est(3,:), 'b-', ...
      x_true_storage(1,:), x_true_storage(2,:), x_true_storage(3,:), 'r--', ...
      y_storage(1,:), y_storage(2,:), y_storage(3,:), 'k.');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
legend('Estimated', 'True', 'Measurements');
title('3D Trajectory');

% Plot curvature, torsion, and tangential velocity estimates
figure('Name', 'Curvature, Torsion, and Velocity Estimation');
subplot(3,1,1);
plot(t, gamma_est, 'b-', t, gamma_true*ones(size(t)), 'r--');
ylabel('Curvature');
legend('Estimated', 'True');
grid on;

subplot(3,1,2);
plot(t, tau_est, 'b-', t, tau_true*ones(size(t)), 'r--');
ylabel('Torsion');
grid on;

subplot(3,1,3);
plot(t, u_est, 'b-', t, u_true*ones(size(t)), 'r--');
ylabel('Tangential Velocity');
xlabel('Time (s)');
grid on;

% Add error analysis
figure('Name', 'Estimation Errors');
error_x = x_est - x_true_storage;
error_y = y_storage - x_true_storage;
subplot(3,1,1);
plot(t, error_y(1,:), 'r-'); hold on
plot(t, error_x(1,:), 'b-'); hold off
ylabel('x error');
legend('Measure Error', 'Estimated Error');
grid on;

subplot(3,1,2);
plot(t, error_y(1,:), 'r-'); hold on
plot(t, error_x(2,:), 'b-'); hold off
ylabel('y error');
grid on;

subplot(3,1,3);
plot(t, error_y(1,:), 'r-'); hold on
plot(t, error_x(3,:), 'b-'); hold off
ylabel('z error');
xlabel('Time (s)');
grid on; 