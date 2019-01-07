x = traj_m;

% kmlwriteline('./plot/traj.kml', x(1:50:end,2),x(1:50:end,3),x(1:50:end,4), 'Color','g', 'Width',2);

vel = zeros(size(x,1),1);
for k=1:size(x,1)
    vel(k) = norm(x(k,5:7));
end

%-------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(x(:,1), x(:,2), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\itL\rm(\circ)')
title('Position')
grid on

subplot(3,1,2)
plot(x(:,1), x(:,3), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\lambda\rm(\circ)')
grid on

subplot(3,1,3)
plot(x(:,1), x(:,4), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\ith\rm(m)')
grid on

%-------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(x(:,1), x(:,5), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\itv_n\rm(m/s)')
title('Velocity')
grid on

subplot(3,1,2)
plot(x(:,1), x(:,6), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\itv_e\rm(m/s)')
grid on

subplot(3,1,3)
plot(x(:,1), x(:,7), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\itv_d\rm(m/s)')
grid on

%-------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(x(:,1), x(:,8), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\psi(\circ)')
title('Attitude')
grid on

subplot(3,1,2)
plot(x(:,1), x(:,9), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\theta(\circ)')
grid on

subplot(3,1,3)
plot(x(:,1), x(:,10), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\gamma(\circ)')
grid on

clearvars x