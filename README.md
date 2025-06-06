# Frenet-Serret and Bishop Invariant Extended Kalman Filters

This repository contains implementations of both Frenet-Serret and Bishop Invariant Extended Kalman Filters (IEKF) for tracking targets in 3D space. The implementations are based on the papers "Tracking the Frenet-Serret frame associated to a highly maneuvering target in 3D" by Pilté et al. and "An Extension to the Frenet-Serret and Bishop Invariant Extended Kalman Filters for Tracking Accelerating Targets" by Gibbs et al.

## Overview

The repository provides two alternative approaches for tracking targets in 3D space:

1. **Frenet-Serret Frame**: Uses curvature and torsion to describe the target's motion
2. **Bishop Frame**: Uses two signed curvatures and includes acceleration in the state vector

Both implementations use the Invariant Extended Kalman Filter (IEKF) framework, which provides better convergence properties than traditional EKFs.

## Frenet-Serret Implementation

The Frenet-Serret implementation models the target's motion using the classical Frenet-Serret formulas, which describe the evolution of the tangent, normal, and binormal vectors along a curve.

### State Vector (Frenet-Serret)
- Position (3x1)
- Rotation matrix (3x3)
- Curvature (scalar)
- Torsion (scalar)
- Tangential velocity (scalar)

### Key Features (Frenet-Serret)
- Handles highly maneuvering targets
- Accounts for unobservable rotations around the velocity vector
- Assumes nearly constant tangential velocity, curvature, and torsion

## Bishop Frame Implementation

The Bishop frame provides an alternative to the Frenet-Serret frame, using two signed curvatures instead of curvature and torsion. This implementation extends the Bishop frame to include acceleration in the state vector.

### State Vector (Bishop)
- Position (3x1)
- Rotation matrix (3x3)
- First Bishop curvature (scalar)
- Second Bishop curvature (scalar)
- Tangential velocity (scalar)
- Tangential acceleration (scalar)

### Key Features (Bishop)
- Better handling of accelerating targets
- More natural representation of certain trajectories
- Improved convergence for certain motion patterns

## Implementation

The repository contains the following main files:

1. Frenet-Serret Implementation:
   - `frenet_serret_ekf.m`: Main filter class
   - `test_frenet_serret_ekf.m`: Test script

2. Bishop Frame Implementation:
   - `bishop_ekf.m`: Main filter class
   - `test_bishop_ekf.m`: Test script

## Usage

### Frenet-Serret Filter
```matlab
% Create filter instance
dt = 0.1;  % Time step
Q_x = eye(3) * 0.1;  % Position process noise
Q_omega = eye(3) * 0.01;  % Angular velocity process noise
Q_gamma = 0.001;  % Curvature process noise
Q_tau = 0.001;  % Torsion process noise
Q_u = 0.1;  % Tangential velocity process noise
R = eye(3) * 0.1;  % Measurement noise

filter = frenet_serret_ekf(dt, Q_x, Q_omega, Q_gamma, Q_tau, Q_u, R);
```
## Output
- Trajectory plots: Show the true target, noisy measurements, and the filter's estimate.
![Screenshot 2025-06-06 185447](https://github.com/user-attachments/assets/ea5140af-71a2-409d-b444-c9f4b34763e2)


### Bishop Filter
```matlab
% Create filter instance
dt = 0.05;  % Time step
Q_x = diag([0.1; 0.1; .2]);     % Position process noise
Q_omega = 0.1 * eye(3);         % Angular velocity process noise
Q_kappa1 = .1;                  % First curvature process noise
Q_kappa2 = .1;                  % Second curvature process noise
Q_u = .1;                       % Tangential velocity process noise
Q_a = .1;                       % Tangential acceleration process noise
R = diag([0.01; 0.01; .02]);    % Measurement noise

filter = bishop_ekf(dt, Q_x, Q_omega, Q_kappa1, Q_kappa2, Q_u, Q_a, R);
```

## Mathematical Background

### Frenet-Serret Formulas
```
d/dt x_t = u T
d/dt T = uκ N
d/dt N = u(-κ T + τ̃ B)
d/dt B = -u τ̃ N
```

### Bishop Frame Formulas
```
d/dt x_t = u T
d/dt T = u(κ₁ M₁ + κ₂ M₂)
d/dt M₁ = -u κ₁ T
d/dt M₂ = -u κ₂ T
```

## References

1. Pilté, M., Bonnabel, S., & Barbaresco, F. (2018). "Tracking the Frenet-Serret frame associated to a highly maneuvering target in 3D". In Proceedings of the IEEE Conference on Decision and Control (CDC).

2. Gibbs, J., Anderson, D., MacDonald, M., & Russell, J. (2022). "An Extension to the Frenet-Serret and Bishop Invariant Extended Kalman Filters for Tracking Accelerating Targets." In Proc. 2022 Sensor Signal Processing for Defence (SSPD), Edinburgh, UK.

3. Bishop, R. L. (1975). "There is more than one way to frame a curve." The American Mathematical Monthly, 82(3), 246-251.

4. Barrau, A., & Bonnabel, S. (2018). "The invariant extended Kalman filter as a stable observer." IEEE Transactions on Automatic Control, 62(4), 1797-1812.

## Dependencies

- MATLAB R2019b or later
- No additional toolboxes required

## License

This project is licensed under the MIT License - see the LICENSE file for details. 
