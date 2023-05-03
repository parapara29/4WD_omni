% Dynamic Model
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
eta0 = [0;0;0]; % actual initial position of the wmr
eta(:,1) = eta0;

% ODE solution (Euler Method)
for i = 1:length(t)
    % Desired Trajectory of the wmr
    eta_des(:,i) = [2-2*cos(0.1*t(i)); 2*sin(0.1*t(i)); 0.1*t(i)];
    eta_des_dot(:,i) = [0.2*sin(0.1*t(i));0.2*cos(0.1*t(i)); 0.1];
  
    x(i) = eta(1,i);
    y(i) = eta(2,i);
    psi(i) = eta(3,i); % generalized/relative orientation of the vehicle
    
    % Jacobain matrix
    J_psi = [cos(psi(i)),-sin(psi(i)), 0; sin(psi(i)), cos(psi(i)), 0; 0, 0, 1];

    % Wheel configuration matrix
    W = [a/2,0,a/2,0; 0,a/2,0,a/2; -(a/(2*l)), -(a/(2*l)),(a/(2*l)),(a/(2*l))];    
    
    % wheel angular velocities
    zeta(:,i) = J_psi'*eta_des_dot(:,i);
    omega(:,i) = pinv(W)*zeta(:,i);
    zeta(:,i) = W*omega(:,i);
     
    eta_dot(:,i) = J_psi * zeta(:,i);

    eta(:,i+1) = eta(:,i) + dt*eta_dot(:,i);
end

% Animation
figure 
plus_v = [-t_edge/2, t_edge/2, t_edge/2, l/2, l/2, t_edge/2, t_edge/2,-t_edge/2, -t_edge/2, -l/2,-l/2, -t_edge/2;
           l/2, l/2, l/2-t_edge/2, l/2-t_edge/2, -(l/2-t_edge/2),-(l/2-t_edge/2), -l/2, -l/2, -(l/2-t_edge/2), -(l/2-t_edge/2),l/2-t_edge/2, l/2-t_edge/2]; % vehicle base
box_w = [-a, a, a, -a, -a; -b/2, -b/2, b/2, b/2, -b/2]; % wheels

for i = 1:length(t)
    R_psi = [cos(psi(i)),-sin(psi(i));sin(psi(i)),+cos(psi(i))];
    
    veh_ani = R_psi * plus_v;
    
    % Wheel positions and orientations
    wheel_1 = R_psi * (R_w1*box_w + [-t_edge/2+a; l/2+b/2]);
    wheel_2 = R_psi * (R_w2*box_w + [l/2+b/2; t_edge/2-a]);
    wheel_3 = R_psi * (R_w3*box_w + [t_edge/2-a; -l/2-b/2]);
    wheel_4 = R_psi * (R_w4*box_w + [-l/2-b/2; -t_edge/2+a]);
    
    fill(veh_ani(1,:)+x(i), veh_ani(2,:)+y(i), 'y');
    hold on
    fill(wheel_1(1,:)+x(i),wheel_1(2,:)+y(i),'r');
    fill(wheel_2(1,:)+x(i),wheel_2(2,:)+y(i),'r');
    fill(wheel_3(1,:)+x(i),wheel_3(2,:)+y(i),'r');
    fill(wheel_4(1,:)+x(i),wheel_4(2,:)+y(i),'r');
    
    plot(x(1:i),y(1:i),'b--', eta_des(1,1:i), eta_des(2,1:i), 'g');
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
legend('x[m]','y[m]','psi[rad]');
xlabel('t[s]');
ylabel('eta[units]');
grid on

figure
plot(t,zeta(1,1:i), 'r-', t,zeta(2,1:i), 'b--', t,zeta(3,1:i), 'g');
legend('u[m/s]','v[m/s]','thetad,[rad/s]');
xlabel('t[s]');
ylabel('zeta[units]');
grid on

figure
plot(t,eta_des(1,1:i)-eta(1,1:i), 'r-', t,eta_des(2,1:i)-eta(2,1:i), 'b--', t,eta_des(3,1:i)-eta(3,1:i), 'g');
legend('xerr[m]','yerr[m]','psierr[rad]');
xlabel('t[s]');
ylabel('error[units]');
grid on

