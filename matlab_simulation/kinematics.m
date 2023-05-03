% Kinematic Model (Land-based mobile robot)
clc; close all;

% Simulation parameters 
dt = 0.1; % step size
ts = 20; % total simulation time
t = 0:dt:ts; % span 

% Physical parameters
l = 0.4; % length of the vehicle (400mm)
t_edge = 0.2; % thickness of the edges (200mm)
a = 0.1; % wheel radius (100mm)
b = 0.05; % wheel thickness (50mm)

% Initial conditions
eta0 = [0; 0; 0]; % actual initial position of the wmr
eta(:,1) = eta0;



% ODE solution (Euler Method)
for i = 1:length(t)
    x(i) = eta(1,i);
    y(i) = eta(2,i);
    psi(i) = eta(3,i); % generalized/relative orientation of the vehicle

     % Update wheel angular velocities for square trajectory
    if i <= length(t)/4
        omega_1 = 1; omega_2 = 0; omega_3 = 1; omega_4 = 0;
    elseif i > length(t)/4 && i <= length(t)/2
        omega_1 = 0; omega_2 = 1; omega_3 = 0; omega_4 = 1;
    elseif i > length(t)/2 && i <= 3*length(t)/4
        omega_1 = -1; omega_2 = 0; omega_3 = -1; omega_4 = 0;
    else
        omega_1 = 0; omega_2 = -1; omega_3 =0; omega_4 =-1;
    end
    % Jacobain matrix
    J_psi = [cos(psi(i)),-sin(psi(i)), 0; sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];

    % Wheel configuration matrix
    W = [a/2,0,a/2,0; 0,a/2,0,a/2; -(a/(2*l)), -(a/(2*l)),(a/(2*l)),(a/(2*l))];
    
    % wheel angular velocities
    zeta(:,i) = W*[omega_1;omega_2;omega_3;omega_4];
     
    eta_dot(:,i) = J_psi * zeta(:,i);

    eta(:,i+1) = eta(:,i) + dt*eta_dot(:,i);
end

% Plus shape coordinates
plus_v = [-t_edge/2, t_edge/2, t_edge/2, l/2, l/2, t_edge/2, t_edge/2,-t_edge/2, -t_edge/2, -l/2,-l/2, -t_edge/2;
           l/2, l/2, l/2-t_edge/2, l/2-t_edge/2, -(l/2-t_edge/2),-(l/2-t_edge/2), -l/2, -l/2, -(l/2-t_edge/2), -(l/2-t_edge/2),l/2-t_edge/2, l/2-t_edge/2];

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

fill(veh_ani(1,:)+x(i),veh_ani(2,:)+y(i),'y');
hold on
fill(wheel_1(1,:)+x(i),wheel_1(2,:)+y(i),'r');
fill(wheel_2(1,:)+x(i),wheel_2(2,:)+y(i),'r');
fill(wheel_3(1,:)+x(i),wheel_3(2,:)+y(i),'r');
fill(wheel_4(1,:)+x(i),wheel_4(2,:)+y(i),'r');

plot(x(1:i),y(1:i),'b-');
set(gca,'fontsize',24)
xlabel('x,[m]');
ylabel('y,[m]');
llim = min(min(x),min(y)) - l;
ulim = max(max(x),max(y)) + l;
axis([llim ulim llim ulim]);
axis square
grid on
pause(0.1)
hold off
end

% Plots
figure
plot(t,eta(1,1:i), 'r-', t,eta(2,1:i), 'b--', t,eta(3,1:i), 'g');
legend('x[m]','y[m]','psi,[rad]');
xlabel('t[s]');
ylabel('eta[units]');
grid on

figure
plot(t,zeta(1,1:i), 'r-', t,zeta(2,1:i), 'b--', t,zeta(3,1:i), 'g');
legend('u[m/s]','v[m]','thetad,[rad/s]');
xlabel('t[s]');
ylabel('eta[units]');
grid on

figure
plot(eta(1,1:i), eta(2,1:i), 'r-')
xlabel('x[m]')
ylabel('y[m]')
grid on
