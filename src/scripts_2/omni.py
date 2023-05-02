import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation

# Simulation parameters
dt = 0.1  # Step size
ts = 5 # Simulation time
t = np.arange(0, ts, dt)  # time span

# Initial conditions
eta0 = np.array([0, 0, 0])  # Initial position and orientation of the vehicle
zeta0 = np.array([0, 0, 0])  # Initial vector of input commands

eta = np.zeros((3, len(t)))
zeta = np.zeros((3, len(t)))

eta[:, 0] = eta0
zeta[:, 0] = zeta0

# Robot parameters
m = 10  # mass of the vehicle
Iz = 0.1  # Inertia of the vehicle

xbc = 0
ybc = 0  # coordinates of mass center

# Wheel configuration parameters
l = 1.0
w = 1.0
thickness = 0.5
wheel_radius = 0.1

phi1 = np.pi / 2
phi2 = 0
phi3 = 3 * np.pi / 2
phi4 = np.pi

# State propagation
for i in range(len(t)):
    u = zeta[0, i]
    v = zeta[1, i]
    r = zeta[2, i]

    # Inertia matrix
    D = np.array([[m, 0, -ybc * m],
                  [0, m, xbc * m],
                  [-ybc * m, xbc * m, Iz + m * (xbc ** 2 + ybc ** 2)]])

    # Other vector
    n_v = np.array([-m * r * (v + xbc * r),
                    m * r * (u - ybc * r),
                    m * r * (xbc * u + ybc * v)])

    # Wheel input matrix (Differential wheel drive)
    Gamma = np.array([[np.cos(phi1), np.cos(phi2), np.cos(phi3), np.cos(phi4)],
                      [np.sin(phi1), np.sin(phi2), np.sin(phi3), np.sin(phi4)],
                      [l, l, l, l]])

    # Wheel inputs (traction forces)
    F1 = 1
    F2 = -0.5
    F3 = -0.5
    F4 = 1
    kappa = np.array([F1, F2, F3, F4])

    # Input vector
    tau = Gamma @ kappa

    # Jacobian matrix
    psi = eta[2, i]
    J_eta = np.array([[np.cos(psi), -np.sin(psi), 0],
                      [np.sin(psi), np.cos(psi), 0],
                      [0, 0, 1]])

    zeta_dot = np.linalg.inv(D) @ (tau - n_v - 0.2 * zeta[:, i])
    zeta[:, i] = zeta[:, i-1] + dt * zeta_dot

    eta[:, i ] = eta[:, i-1] + dt * J_eta @ (zeta[:, i] + dt / 2 * zeta_dot)

# Plotting and Animation
def plot_robot(frame):
    i = frame
    ax.clear()
    psi = eta[2, i]

    # Robot body coordinates
    robot_body = np.array([[-l / 2, l / 2, l / 2, -l / 2, -l / 2, -l / 2, 0, 0, l / 2, l / 2],
                           [-w / 2, -w / 2, w / 2, w / 2, -w / 2, w / 2, w / 2, -w / 2, -w / 2, w / 2]])

    # Wheel coordinates
    wheel = 0.5 * np.array([[-thickness, thickness, thickness, -thickness, -thickness],
                            [-thickness, -thickness, thickness, thickness, -thickness]])

    R_psi = np.array([[np.cos(psi), -np.sin(psi)],
                      [np.sin(psi), np.cos(psi)]])

    wheel_positions = [np.array([l / 2 * np.cos(phi1 - np.pi / 2), l / 2 * np.sin(phi1 - np.pi / 2)]),
                       np.array([l / 2 * np.cos(phi2 - np.pi / 2), l / 2 * np.sin(phi2 - np.pi / 2)]),
                       np.array([l / 2 * np.cos(phi3 - np.pi / 2), l / 2 * np.sin(phi3 - np.pi / 2)]),
                       np.array([l / 2 * np.cos(phi4 - np.pi / 2), l / 2 * np.sin(phi4 - np.pi / 2)])]

    robot_pos = R_psi @ robot_body
    ax.add_patch(patches.Polygon(robot_pos.T + eta[:2, i], closed=True, color='g'))

    for j, phi in enumerate([phi1, phi2, phi3, phi4]):
        R_phi = np.array([[np.cos(phi), -np.sin(phi)],
                          [np.sin(phi), np.cos(phi)]])
        wheel_pos = R_psi @ (R_phi @ wheel + wheel_positions[j].reshape(-1, 1))
        ax.add_patch(patches.Polygon(wheel_pos.T + eta[:2, i], closed=True, color='r'))

    l_lim = np.min(eta[:2, :])
    u_lim = np.max(eta[:2, :])
    ax.set_xlim(-0.5 + l_lim, 0.5 + u_lim)
    ax.set_ylim(-0.5 + l_lim, 0.5 + u_lim)
    ax.set_aspect('equal', adjustable='box')

    ax.plot(eta[0, :i + 1], eta[1, :i + 1], 'b-')
    ax.grid(True)
    ax.set_xlabel('x, [m]')
    ax.set_ylabel('y, [m]')


fig, ax = plt.subplots()
ani = FuncAnimation(fig, plot_robot, frames=len(t), interval=100, repeat=False)
plt.show()

# Plotting functions
plt.figure()
plt.plot(t, eta[0, :], 'r-', t, eta[1, :], 'b-.', t, eta[2, :], 'k--', linewidth=2)
plt.legend(['x, [m]', 'y, [m]', r'$\psi$, [rad]'])
plt.grid(True)
plt.xlabel('t, [s]')
plt.ylabel('eta, [units]')
plt.show()

plt.figure()
plt.plot(t, zeta[0, :], 'r-', t, zeta[1, :], 'b-.', t, zeta[2, :], 'k--', linewidth=2)
plt.legend(['u, [m/s]', 'v, [m/s]', 'r, [rad/s]'])
plt.grid(True)
plt.xlabel('t, [s]')
plt.ylabel('zeta, [units]')
plt.show()