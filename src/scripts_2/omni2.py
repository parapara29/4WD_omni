import numpy as np
import matplotlib.pyplot as plt

# Simulation parameters
dt = 0.1  # Step size
ts = 5  # Simulation time
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
l = 0.3
phi1 = 90 * np.pi / 180
phi2 = 210 * np.pi / 180
phi3 = 330 * np.pi / 180

# State propagation
for i in range(len(t)):
    u, v, r = zeta[:, i]

    # Inertia matrix
    D = np.array([[m, 0, -ybc * m],
                  [0, m, xbc * m],
                  [-ybc * m, xbc * m, Iz + m * (xbc ** 2 + ybc ** 2)]])

    # Other vector
    n_v = np.array([-m * r * (v + xbc * r),
                    m * r * (u - ybc * r),
                    m * r * (xbc * u + ybc * v)])

    # Wheel input matrix (Differential wheel drive)
    Gamma = np.array([[np.cos(phi1), np.cos(phi2), np.cos(phi3)],
                      [np.sin(phi1), np.sin(phi2), np.sin(phi3)],
                      [l, l, l]])

    # Wheel inputs (traction forces)
    F1 = 1
    F2 = -0.5
    F3 = -0.5
    kappa = np.array([F1, F2, F3])

    # Input vector
    tau = np.dot(Gamma, kappa)

    # Jacobian matrix
    psi = eta[2, i]
    J_eta = np.array([[np.cos(psi), -np.sin(psi), 0],
                      [np.sin(psi), np.cos(psi), 0],
                      [0, 0, 1]])

    zeta_dot = np.linalg.inv(D).dot(tau - n_v - 0.2 * zeta[:, i])
    zeta[:, i + 1] = zeta[:, i] + dt * zeta_dot

    eta[:, i + 1] = eta[:, i] + dt * J_eta.dot(zeta[:, i] + dt / 2 * zeta_dot)

# Animation (mobile robot motion animation)
l = 0.6  # length of the mobile robot
w = 0.6  # width of the mobile robot
# Mobile robot coordinates
mr_co = np.array([[-l/2, l/2, l/2, -l/2, -l/2],
                  [-w/2, -w/2, w/2, w/2, -w/2]])

wheel_b = 0.5 * np.array([[-0.2, 0.2, -0.2, -0.2, -0.2],
[-0.05, -0.05, 0.05, 0.05, -0.05]])

fig = plt.figure()
for i in range(len(t)):
    psi = eta[2, i]
    R_psi = np.array([[np.cos(psi), -np.sin(psi)],
    [np.sin(psi), np.cos(psi)]]) # rotation matrix
    R_1 = np.array([[np.cos(phi1), -np.sin(phi1)],
    [np.sin(phi1), np.cos(phi1)]]) # rotation matrix of wheel 1
    w_pos1 = np.dot(R_psi, np.dot(R_1, wheel_b) + np.array([[0.2 * np.cos(phi1 - np.pi / 2)],
    [0.2 * np.sin(phi1 - np.pi / 2)]]))
    R_2 = np.array([[np.cos(phi2), -np.sin(phi2)],
    [np.sin(phi2), np.cos(phi2)]]) # rotation matrix of wheel 2
    w_pos2 = np.dot(R_psi, np.dot(R_2, wheel_b) + np.array([[0.2 * np.cos(phi2 - np.pi / 2)],
    [0.2 * np.sin(phi2 - np.pi / 2)]]))
    R_3 = np.array([[np.cos(phi3), -np.sin(phi3)],
    [np.sin(phi3), np.cos(phi3)]]) # rotation matrix of wheel 3
    w_pos3 = np.dot(R_psi, np.dot(R_3, wheel_b) + np.array([[0.2 * np.cos(phi3 - np.pi / 2)],
    [0.2 * np.sin(phi3 - np.pi / 2)]]))
    v_pos = np.dot(R_psi, mr_co)
    plt.fill(v_pos[0, :] + eta[0, i], v_pos[1, :] + eta[1, i], 'g')
    plt.fill(w_pos1[0, :] + eta[0, i], w_pos1[1, :] + eta[1, i], 'r')
    plt.fill(w_pos2[0, :] + eta[0, i], w_pos2[1, :] + eta[1, i], 'r')
    plt.fill(w_pos3[0, :] + eta[0, i], w_pos3[1, :] + eta[1, i], 'r')
    l_lim = np.min(np.min(eta[:2, :]))
    u_lim = np.max(np.max(eta[:2, :]))
    plt.axis([-0.5 + l_lim, 0.5 + u_lim, -0.5 + l_lim, 0.5 + u_lim])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.plot(eta[0, :i+1], eta[1, :i+1], 'b-')
    plt.xlabel('x,[m]')
    plt.ylabel('y,[m]')
    plt.grid(True)
    plt.pause(0.1)
    plt.clf()

plt.figure()
plt.plot(t, eta[0, :], 'r-', label='x,[m]')
plt.plot(t, eta[1, :], 'b-.', label='y,[m]')
plt.plot(t, eta[2, :], 'k--', label='\psi,[rad]')
plt.legend()
plt.grid(True)
plt.xlabel('t,[s]')
