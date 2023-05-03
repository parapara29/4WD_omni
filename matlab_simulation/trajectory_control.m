% Dynamic model
clc; close all;

% Simulation Parameters
dt = 0.1; % step size
ts = 32; % simulation time
t = 0:dt:ts; % time span

%Initial conditions
eta0 = [0;0;0]; % actual initial position of the wmr
zeta0 = [0;0;0]; % inital motion of the wmr

eta(:,1) = eta0;
zeta(:,1) = zeta0;

% Physical parameters
w = 0.5; % width of the vehicle
l = 0.5; % length of the vehicle
a = 0.1; % wheel radius
b = 0.05; % wheel width
h = l/4;
d = w/2; % half of the distance between two wheel centers

%Vehicle parameters
m = 60; Iz = 0.5007;
xbc = 0; ybc = 0;

% controller Parameters
Kp = 9;
Kd = 6;

% ODE solution (Euler Method)
for i = 1:length(t)
    x(i) = eta(1,i); y(i) = eta(2,i); psi(i) = eta(3,i);
    u = zeta(1,i); v = zeta(2,i); r = zeta(3,i);

    % Desired Values
    eta_des(:,i) = [2*sin(0.1*t(i));2-2*cos(0.1*t(i));0.1*t(i)];
    eta_des_dot(:,i) = [0.2*cos(0.1*t(i));0.2*sin(0.1*t(i));0.1];
    eta_des_ddot(:,i) = [-0.02*sin(0.1*t(i));0.02*cos(0.1*t(i));0];
    
    % Inertial Part
    D = [m,0,0; 0,m,0; 0,0,Iz];
    n_zeta = [-m*v*r;m*u*r;0];
    f_zeta = 0*0.3*zeta(:,i); % without friction
    J_eta = [cos(psi(i)),-sin(psi(i)),0; sin(psi(i)),+cos(psi(i)),0; 0,0,1];
    
    tau(:,i) = D*J_eta'*(eta_des_ddot(:,i)+Kp*(eta_des(:,i)-eta(:,i)) +Kd*(eta_des_dot(:,i)-J_eta*zeta(:,i))) + n_zeta;
    
    % Wheel configuration matrix
    gamma = [1,0,1,0; 0,1,0,1; -l/2,-l/2,l/2,l/2];
    kappa(:,i) = pinv(gamma)*tau(:,i);

    zeta_dot(:,i) = inv(D)*(gamma*kappa(:,i) - n_zeta - f_zeta);
    zeta(:,i+1) = zeta(:,i) + zeta_dot(:,i)*dt;
    
    eta_dot(:,i) = J_eta*(zeta(:,i)+zeta_dot(:,i)*dt/2);
    eta(:,i+1) = eta(:,i) + eta_dot(:,i)*dt;
end

% Animation
box_v = [-l/2,l/2,l/2,-l/2,-l/2; -w/2,-w/2,w/2,w/2,-w/2]; % vehicle base

box_w = [-a, a, a, -a, -a; -b/2, -b/2, b/2, b/2, -b/2]; % wheels

for i = 1:length(t)
    
    R_psi = [cos(psi(i)),-sin(psi(i));sin(psi(i)),+cos(psi(i))];
    % Configuration 1
    R_w1 = [cosd(0),-sind(0);sind(0),+cosd(0)];
    R_w2 = [cosd(90),-sind(90);sind(90),+cosd(90)];
    R_w3 = [cosd(0),-sind(0);sind(0),+cosd(0)];
    R_w4 = [cosd(90),-sind(90);sind(90),+cosd(90)];

        
    veh_ani = R_psi * box_v;
    % Configuration 1
    wheel_1 = R_psi * (R_w1*box_w+[(l/2)*cosd(0);(l/2)*sind(90)]);
    wheel_2 = R_psi * (R_w2*box_w+[-(l/2)*cosd(0);(l/2)*sind(90)]);
    wheel_3 = R_psi * (R_w3*box_w+[-(l/2)*cosd(0);-(l/2)*sind(90)]);
    wheel_4 = R_psi * (R_w4*box_w+[(l/2)*cosd(0);-(l/2)*sind(90)]);
 
    fill(veh_ani(1,:)+x(i),veh_ani(2,:)+y(i),'y');
    hold on
    fill(wheel_1(1,:)+x(i),wheel_1(2,:)+y(i),'r');
    fill(wheel_2(1,:)+x(i),wheel_2(2,:)+y(i),'r');
    fill(wheel_3(1,:)+x(i),wheel_3(2,:)+y(i),'r');
    fill(wheel_4(1,:)+x(i),wheel_4(2,:)+y(i),'r');
    
    plot(x(1:i),y(1:i),'b--',eta_des(1,1:i),eta_des(2,1:i),'g-');
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
plot(t,eta_des(1,1:i)-eta(1,1:i), 'r-', t,eta_des(2,1:i)-eta(2,1:i), 'b--', t,eta_des(3,1:i)-eta(3,1:i), 'g');
legend('xerr[m]','yerr[m]','psierr[rad]');
grid on
xlabel('t[units]');
ylabel('error[units]');

figure
plot(t,eta(1,1:i), 'r-', t,eta(2,1:i), 'b--', t,eta(3,1:i), 'g');
grid on
legend('x,[m]','y,[m]','\psi,[rad]');
xlabel('t[s]');
ylabel('eta[units]');

