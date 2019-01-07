% 认为速度是常数模型
T = 5;

t = 0:0.1:T;
vehicle = zeros(length(t),10);
for k=1:length(t)
    vehicle(k,:) = vehicle_traj(t(k));
end

target = zeros(length(t),7);
for k=1:length(t)
    target(k,:) = target_traj(t(k));
end

figure
plot3(vehicle(:,2),vehicle(:,3),vehicle(:,4), 'LineWidth',2)
hold on
grid on
axis equal
plot3(target(:,2),target(:,3),target(:,4), 'LineWidth',2)

%%
dt = target(2,1);
X = zeros(6,1);
Phi = eye(6);
Phi(1,2) = dt;
Phi(3,4) = dt;
Phi(5,6) = dt;
Gamma = zeros(6,3);
Gamma(1,1) = 0.5*dt^2;
Gamma(2,1) = dt;
Gamma(3,2) = 0.5*dt^2;
Gamma(4,2) = dt;
Gamma(5,3) = 0.5*dt^2;
Gamma(6,3) = dt;
P = diag([50,30, 50,30, 50,30])^2;
Q = diag([1, 1, 1])^2;
R = diag([1/60/180*pi, 1/60/180*pi, 5])^2;

%%
figure
h1 = plot3(vehicle(:,2),vehicle(:,3),vehicle(:,4), 'Color','b', 'LineWidth',1); %vehicle trajectory
hold on
grid on
axis equal
h2 = plot3(target(:,2),target(:,3),target(:,4), 'Color','r', 'LineWidth',1); %target trajectory
h3 = plot3(vehicle(1,2),vehicle(1,3),vehicle(1,4), 'Color','b', 'Marker','.', 'MarkerSize',16); %vehicle point
h4 = plot3(target(1,2),target(1,3),target(1,4), 'Color','r', 'Marker','.', 'MarkerSize',16); %target point
h5 = plot3([vehicle(1,2),target(1,2)],[vehicle(1,3),target(1,3)],[vehicle(1,4),target(1,4)], 'Color','m', 'LineStyle','--');
h6 = plot3(X(1),X(3),X(5), 'Color','c', 'Marker','.', 'MarkerSize',16); %estimation pint
h7 = plot3(X(1),X(3),X(5), 'Color','c', 'LineWidth',1); %estimation trajectory
h8 = plot3([vehicle(1,2),X(1)],[vehicle(1,3),X(3)],[vehicle(1,4),X(5)], 'Color','g', 'LineStyle','--');

%%
n = T/dt+1;
filter_X = zeros(n,6);
filter_X(1,:) = X';
filter_P = zeros(n,6);
filter_P(1,:) = sqrt(diag(P))';

for k=2:n
    Rv = vehicle(k,2:4);
    Rt = target(k,2:4);
    r = Rt - Rv;
    
    Z = zeros(3,1);
    Z(1) = atan2(r(2),r(1)) + randn(1)*(1/60/180*pi);
    Z(2) = atan2(r(3),norm(r(1:2))) + randn(1)*(1/60/180*pi);
    Z(3) = norm(r) + randn(1)*5;
    
    X = Phi*X;
    P = Phi*P*Phi' + Gamma*Q*Gamma';
        
    Rt = [X(1),X(3),X(5)];
    r = Rt - Rv;
    
    H = zeros(3,6);
    r1 = r(1)^2 + r(2)^2;
    r2 = r(1)^2 + r(2)^2 + r(3)^2;
    H(1,1) = -r(2)/r1;
    H(1,3) =  r(1)/r1;
    H(2,1) = -r(1)*r(3)/r2/sqrt(r1);
    H(2,3) = -r(2)*r(3)/r2/sqrt(r1);
    H(2,5) =   sqrt(r1)/r2;
    H(3,1) = r(1)/sqrt(r2);
    H(3,3) = r(2)/sqrt(r2);
    H(3,5) = r(3)/sqrt(r2);
    
    Zp = zeros(3,1);
    Zp(1) = atan2(r(2),r(1));
    Zp(2) = atan2(r(3),norm(r(1:2)));
    Zp(3) = norm(r);
    
    K = P*H' / (H*P*H'+R);
    X = X + K*(Z-Zp);
    P = (eye(length(X))-K*H)*P;
    P = (P+P')/2;
    
    filter_X(k,:) = X';
    filter_P(k,:) = sqrt(diag(P))';

    %---------------------------------------------------------------------%
    Rv = vehicle(k,2:4);
    Rt = target(k,2:4);
    set(h3, 'XData',Rv(1), 'YData',Rv(2), 'ZData',Rv(3)); %vehicle point
    set(h4, 'XData',Rt(1), 'YData',Rt(2), 'ZData',Rt(3)); %target point
    set(h5, 'XData',[Rv(1),Rt(1)], 'YData',[Rv(2),Rt(2)], 'ZData',[Rv(3),Rt(3)]);
    set(h6, 'XData',X(1), 'YData',X(3), 'ZData',X(5)); %estimation pint
    set(h7, 'XData',filter_X(1:k,1), 'YData',filter_X(1:k,3), 'ZData',filter_X(1:k,5)); %estimation trajectory
    set(h8, 'XData',[Rv(1),X(1)], 'YData',[Rv(2),X(3)], 'ZData',[Rv(3),X(5)]);
    pause(dt);
end

%%
t = 0:dt:T;

figure
subplot(3,1,1)
plot(t, filter_X(:,1)-target(:,2), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltax\rm(m)')
title('目标位置估计误差')

subplot(3,1,2)
plot(t, filter_X(:,3)-target(:,3), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltay\rm(m)')

subplot(3,1,3)
plot(t, filter_X(:,5)-target(:,4), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltaz\rm(m)')

figure
subplot(3,1,1)
plot(t, filter_X(:,2), 'LineWidth',2)
grid on
hold on
plot(t, target(:,5), 'LineStyle','--')
xlabel('\itt\rm(s)')
ylabel('\itv_x\rm(m/s)')
title('目标速度估计')

subplot(3,1,2)
plot(t, filter_X(:,4), 'LineWidth',2)
grid on
hold on
plot(t, target(:,6), 'LineStyle','--')
xlabel('\itt\rm(s)')
ylabel('\itv_y\rm(m/s)')

subplot(3,1,3)
plot(t, filter_X(:,6), 'LineWidth',2)
grid on
hold on
plot(t, target(:,7), 'LineStyle','--')
xlabel('\itt\rm(s)')
ylabel('\itv_z\rm(m/s)')

figure
subplot(3,1,1)
plot(t, filter_X(:,2)-target(:,5), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltav_x\rm(m/s)')
title('目标速度估计误差')

subplot(3,1,2)
plot(t, filter_X(:,4)-target(:,6), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltav_y\rm(m/s)')

subplot(3,1,3)
plot(t, filter_X(:,6)-target(:,7), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltav_z\rm(m/s)')

%%
function traj = vehicle_traj(t)
    p0 = [-2500, -2500, 5125];
    v0 = [500, 500, -1000]*0;
    a = [0, 0, -10]*0;
    traj = zeros(1,10);
    traj(1) = t;
    traj(2:4) = p0 + v0*t + 0.5*a*t^2;
    traj(5:7) = v0 + a*t;
    traj(8:10) = a;
end

function traj = target_traj(t)
    p0 = [100, 100, 100];
    v0 = [100, 50, 30];
    a = [0, 0, 0];
    traj = zeros(1,7);
    traj(1) = t;
    traj(2:4) = p0 + v0*t + 0.5*a*t^2;
    traj(5:7) = v0 + a*t;
end

% function traj = target_traj(t)
%     p0 = [100, 100, 0];
%     v = 100;
%     r = 500;
%     traj = zeros(1,7);
%     traj(1) = t;
%     traj(2) = p0(1) + r*cos(v/r*t) - r;
%     traj(3) = p0(2) + r*sin(v/r*t);
%     traj(4) = p0(3);
%     traj(5) = -v*sin(v/r*t);
%     traj(6) =  v*cos(v/r*t);
%     traj(7) = 0;
% end