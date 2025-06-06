% Test script for Bishop EKF implementation

% Parameters
dt = 0.05;  % Time step
T = 10;     % Simulation time
N = T/dt;   % Number of steps

% Process noise parameters
Q_x = diag([0.1; 0.1; .2]);     % Position process noise
Q_omega = 0.1 * eye(3);            % Angular velocity process noise
Q_kappa1 = .1;                     % First curvature process noise
Q_kappa2 = .1;                     % Second curvature process noise
Q_u = .1;                          % Tangential velocity process noise
Q_a = .1;                          % Tangential acceleration process noise

% Measurement noise
R = diag([0.01; 0.01; .02]);          % Measurement noise

% Initialize EKF
ekf = bishop_ekf(dt, Q_x, Q_omega, Q_kappa1, Q_kappa2, Q_u, Q_a, R);

% Initialize true state
x_true = zeros(3,1);
R_true = eye(3);
kappa1_true = 0.1;  % Constant first curvature
kappa2_true = 0.05; % Constant second curvature
u_true = 10;        % Initial tangential velocity
a_true = 0.01;       % Constant tangential acceleration

% Storage for plotting
x_est = zeros(3,N);
x_true_storage = zeros(3,N);
y_storage = zeros(3,N);     % Store measurements
kappa1_est = zeros(1,N);    % Store first curvature estimates
kappa2_est = zeros(1,N);    % Store second curvature estimates
u_est = zeros(1,N);         % Store tangential velocity estimates
a_est = zeros(1,N);         % Store tangential acceleration estimates
t = 0:dt:T-dt;

% Main loop
for k = 1:N
    % True state propagation
    omega = [0; -kappa2_true; kappa1_true];  % Bishop Darboux vector
    Omega = skew(omega);
    R_true = R_true * expm(Omega * dt);
    
    % Update velocity with acceleration
    u_true = u_true + a_true * dt;
    
    % Update position
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
    kappa1_est(k) = ekf.kappa1_hat;
    kappa2_est(k) = ekf.kappa2_hat;
    u_est(k) = ekf.u_hat;
    a_est(k) = ekf.a_hat;
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

% Plot curvature, velocity and acceleration estimates
figure('Name', 'State Estimation');
subplot(4,1,1);
plot(t, kappa1_est, 'b-', t, kappa1_true*ones(size(t)), 'r--');
ylabel('First Curvature');
legend('Estimated', 'True');
grid on;

subplot(4,1,2);
plot(t, kappa2_est, 'b-', t, kappa2_true*ones(size(t)), 'r--');
ylabel('Second Curvature');
grid on;

subplot(4,1,3);
plot(t, u_est, 'b-', t, u_true + a_true*t, 'r--');
ylabel('Tangential Velocity');
grid on;

subplot(4,1,4);
plot(t, a_est, 'b-', t, a_true*ones(size(t)), 'r--');
ylabel('Tangential Acceleration');
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
plot(t, error_y(2,:), 'r-'); hold on
plot(t, error_x(2,:), 'b-'); hold off
ylabel('y error');
grid on;

subplot(3,1,3);
plot(t, error_y(3,:), 'r-'); hold on
plot(t, error_x(3,:), 'b-'); hold off
ylabel('z error');
xlabel('Time (s)');
grid on; 