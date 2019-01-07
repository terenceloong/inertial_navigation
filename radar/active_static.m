T = 5;

t = 0:0.1:T;
vehicle = zeros(length(t),10);
for k=1:length(t)
    vehicle(k,:) = vehicle_traj(t(k));
end

target = [0, 0, 0];

figure
plot3(vehicle(:,2),vehicle(:,3),vehicle(:,4), 'LineWidth',2)
hold on
grid on
axis equal
plot3(target(1),target(2),target(3), 'Marker','.', 'MarkerSize',16)

%%
dt = vehicle(2,1);
X = [vehicle(1,2);vehicle(1,5);vehicle(1,3);vehicle(1,6);vehicle(1,4);vehicle(1,7)];
dX = [400;2; 400;-2; 400;1.5];
X = X + dX;
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
P = diag([500,1, 500,1, 500,1])^2;
Q = diag([0.01, 0.01, 0.01])^2;
R = diag([1/60/180*pi, 1/60/180*pi, 5])^2;

%%
figure
h1 = plot3(vehicle(:,2),vehicle(:,3),vehicle(:,4), 'Color','b', 'LineWidth',1); %vehicle trajectory
hold on
grid on
axis equal
h2 = plot3(target(1),target(2),target(3), 'Color','r', 'Marker','.', 'MarkerSize',16); %target point
h3 = plot3(vehicle(1,2),vehicle(1,3),vehicle(1,4), 'Color','b', 'Marker','.', 'MarkerSize',16); %vehicle point
h4 = plot3([vehicle(1,2),target(1)],[vehicle(1,3),target(2)],[vehicle(1,4),target(3)], 'Color','m', 'LineStyle','--');
h5 = plot3(X(1),X(3),X(5), 'Color','c', 'Marker','.', 'MarkerSize',16); %navigation point
h6 = plot3(X(1),X(3),X(5), 'Color','c', 'LineWidth',1); %navigation trajectory
h7 = plot3([X(1),target(1)],[X(3),target(2)],[X(5),target(3)], 'Color','g', 'LineStyle','--');

%%
n = T/dt+1;
filter_X = zeros(n,6);
filter_X(1,:) = X';
filter_P = zeros(n,6);
filter_P(1,:) = sqrt(diag(P))';

for k=2:n
    acc = vehicle(k,8:10)' + randn(3,1)*0.01;
    
    X = Phi*X + Gamma*acc;
    P = Phi*P*Phi' + Gamma*Q*Gamma';
    
%     if mod(k-1,10)==0
        r = vehicle(k,2:4) - target;

        Z = zeros(3,1);
        Z(1) = atan2(r(2),r(1)) + randn(1)*(1/60/180*pi);
        Z(2) = atan2(r(3),norm(r(1:2))) + randn(1)*(1/60/180*pi);
        Z(3) = norm(r) + randn(1)*5;
        
        r = [X(1),X(3),X(5)] - target;

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
%     end
    
    filter_X(k,:) = X';
    filter_P(k,:) = sqrt(diag(P))';

    %---------------------------------------------------------------------%
    set(h3, 'XData',vehicle(k,2), 'YData',vehicle(k,3), 'ZData',vehicle(k,4)); %vehicle point
    set(h4, 'XData',[vehicle(k,2),target(1)], 'YData',[vehicle(k,3),target(2)], 'ZData',[vehicle(k,4),target(3)]);
    set(h5, 'XData',X(1), 'YData',X(3), 'ZData',X(5)); %navigation point
    set(h6, 'XData',filter_X(1:k,1), 'YData',filter_X(1:k,3), 'ZData',filter_X(1:k,5)); %navigation trajectory
    set(h7, 'XData',[X(1),target(1)], 'YData',[X(3),target(2)], 'ZData',[X(5),target(3)]);
    pause(dt);
end

%%
t = 0:dt:T;

figure %position error
subplot(3,1,1)
plot(t, filter_X(:,1)-vehicle(:,2), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltax\rm(m)')
title('导航位置误差')

subplot(3,1,2)
plot(t, filter_X(:,3)-vehicle(:,3), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltay\rm(m)')

subplot(3,1,3)
plot(t, filter_X(:,5)-vehicle(:,4), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltaz\rm(m)')

figure %velocity error
subplot(3,1,1)
plot(t, filter_X(:,2)-vehicle(:,5), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltav_x\rm(m/s)')
title('导航速度误差')

subplot(3,1,2)
plot(t, filter_X(:,4)-vehicle(:,6), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltav_y\rm(m/s)')

subplot(3,1,3)
plot(t, filter_X(:,6)-vehicle(:,7), 'LineWidth',2)
grid on
xlabel('\itt\rm(s)')
ylabel('\it\deltav_z\rm(m/s)')

%%
function traj = vehicle_traj(t)
    p0 = [-2500, -2500, 5125] + [-1000, -1000, 500];
    v0 = [500, 500, -1000];
    a = [0, 0, -10];
%     p0 = [-2500, -2500, 800];
%     v0 = [50, 50, 0];
%     a = [0, 0, 0];
    traj = zeros(1,10);
    traj(1) = t;
    traj(2:4) = p0 + v0*t + 0.5*a*t^2;
    traj(5:7) = v0 + a*t;
    traj(8:10) = a;
end