% Dynamic Model

clc; close all;
% Simulation parameters
ts = 20; % total simulation time
tspan = [0 ts]; % span

% Physical parameters (updated)
l = 0.4; % length of the vehicle (400mm)
t_edge = 0.2; % thickness of the edges (200mm)
a = 0.1; % wheel radius (100mm)
b = 0.05; % wheel thickness (50mm)

% Initial conditions
eta0 = [0;0;0]; % actual initial position of the wmr
zeta0 = [0;0;0]; % inital motion of the wmr
initial_conditions = [eta0; zeta0];

% Mwr Physical Parameters
m = 10; % maximum mass of wmr
Iz = 0.9;

xbc = 0;
ybc = 0;

% ODE solution (ODE45)
[t, state] = ode45(@(t, state) wmr_dynamics(t, state, l, a, b, m, Iz, xbc, ybc), tspan, initial_conditions);

% Extract state variables
eta = state(:, 1:3)';
zeta = state(:, 4:6)';

% Desired trajectory
%eta_des = [0.1*t.^2;0.1*t.^2;0];
%eta_des_dot = [0.2*t;0.2*t;0];
%eta_des_ddot = [0.2;0.2;0];
eta_des = [2*sin(0.1*t);2-2*cos(0.1*t);0.1*t];
eta_des_dot = [0.2*cos(0.1*t);0.2*sin(0.1*t);0.1];
eta_des_ddot = [-0.02*sin(0.1*t);0.02*cos(0.1*t);0];

% Animation
% Plus shape coordinates (updated)
plus_v = [-t_edge/2, t_edge/2, t_edge/2, l/2, l/2, t_edge/2, t_edge/2,-t_edge/2, -t_edge/2, -l/2,-l/2, -t_edge/2;
           l/2, l/2, l/2-t_edge/2, l/2-t_edge/2, -(l/2-t_edge/2),-(l/2-t_edge/2), -l/2, -l/2, -(l/2-t_edge/2), -(l/2-t_edge/2), l/2-t_edge/2, l/2-t_edge/2]; % vehicle shape (updated)
box_w = [-a, a, a, -a, -a; -b/2, -b/2, b/2, b/2, -b/2]; % wheels
% Animation
for i = 1:length(t)
    R_psi = [cos(psi(i)),-sin(psi(i));sin(psi(i)),+cos(psi(i))];
    R_w1 = [cosd(0),-sind(0);sind(0),+cosd(0)];
    R_w2 = [cosd(90),-sind(90);sind(90),+cosd(90)];
    R_w3 = [cosd(0),-sind(0);sind(0),+cosd(0)];
    R_w4 = [cosd(90),-sind(90);sind(90),+cosd(90)];
    veh_ani = R_psi * plus_v;

    % Wheel positions and orientations
    wheel_1 = R_psi * (R_w1*box_w + [-t_edge/2+a; l/2+b/2]);
    wheel_2 = R_psi * (R_w2*box_w + [l/2+b/2; t_edge/2-a]);
    wheel_3 = R_psi * (R_w3*box_w + [t_edge/2-a; -l/2-b/2]);
    wheel_4 = R_psi * (R_w4*box_w + [-l/2-b/2; -t_edge/2+a]);
    
    fill(veh_ani(1,:)+eta(1,i),veh_ani(2,:)+eta(2,i),'y');
    hold on
    fill(wheel_1(1,:)+eta(1,i),wheel_1(2,:)+eta(2,i),'r');
    fill(wheel_2(1,:)+eta(1,i),wheel_2(2,:)+eta(2,i),'r');
    fill(wheel_3(1,:)+eta(1,i),wheel_3(2,:)+eta(2,i),'r');
    fill(wheel_4(1,:)+eta(1,i),wheel_4(2,:)+eta(2,i),'r');
   plot(eta(1,1:end),eta(2,1:end),'b--',eta_des(1,1:end), eta_des(2,1:end), 'g');

    set(gca,'fontsize',24)
    xlabel('x,[m]');
    ylabel('y,[m]');
    llim = min(min(eta(1,:)),min(eta(2,:))) - l;
    ulim = max(max(eta(1,:)),max(eta(2,:))) + l;
    axis([llim ulim llim ulim]);
    axis square
    grid on
    pause(0.1)
    hold off
end

% Calculate errors
x_error = eta_des(1,:) - eta(1,1:length(t))';
y_error = eta_des(2,:) - eta(2,1:length(t))';
psi_error = eta_des(3,:) - eta(3,1:length(t))';

% Plots
figure
plot(t,eta(1,:), 'r-', t,eta(2,:), 'b--', t,eta(3,:), 'g');
legend('x[m]','y[m]','psi[rad]');
xlabel('t[s]');
ylabel('eta[units]');
grid on

figure
plot(t, x_error, 'r-', t, y_error, 'b--', t, psi_error, 'g');
legend('xerr[m]', 'yerr[m]', 'psierr[rad]');
xlabel('t[s]');
ylabel('error[units]');
grid on


function [dstatedt] = wmr_dynamics(t, state, l, a, b, m, Iz, xbc, ybc)
    eta = state(1:3);
    zeta = state(4:6);
    
    u = zeta(1); v = zeta(2); r = zeta(3);
    
    % Inertia Matrix
    D = [m, 0, 0; 0, m, 0; 0, 0, Iz];
    
    % Other Forces
    n_v = [-m*r*(v+xbc*r); m*r*(u-ybc*r); m*r*(xbc*u+ybc*v)];
    
    % Control law (PD control)
    Kp = 100 * eye(3);
    Kd = 10 * eye(3);
    eta_des = [2*sin(0.1*t); 2-2*cos(0.1*t); 0.1*t];
    eta_des_dot = [0.2*cos(0.1*t); 0.2*sin(0.1*t); 0.1];
    eta_des_ddot = [-0.02*sin(0.1*t); 0.02*cos(0.1*t); 0];
    
    eta_err = eta_des - eta;
    eta_err_dot = eta_des_dot - zeta;
    
    tau = D*eta_des_ddot + n_v + Kd*eta_err_dot + Kp*eta_err;
    
    % Wheel configuration matrix(forces and moments)
    gamma = [1, 1, 1, 1; 1, 1, 1, 1; l/2, -l/2, l/2, -l/2];
    
    % Input forces
    kappa = pinv(gamma) * tau;
    
    % Jacobian matrix
    J_psi = [cos(eta(3)), -sin(eta(3)), 0; sin(eta(3)), cos(eta(3)), 0; 0, 0, 1];
    
    % State derivative
    deta_dt = J_psi * zeta;
    dzeta_dt = pinv(D) * (gamma * kappa - n_v);
    
    dstatedt = [deta_dt; dzeta_dt];
end
