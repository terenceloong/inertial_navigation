x = nav_error;

%-------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(x(:,1), x(:,2), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\itL\rm(m)')
title('Position error')
grid on

subplot(3,1,2)
plot(x(:,1), x(:,3), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\lambda(m)')
grid on

subplot(3,1,3)
plot(x(:,1), x(:,4), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\ith\rm(m)')
grid on

%-------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(x(:,1), x(:,5), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_n\rm(m/s)')
title('Velocity error')
grid on

subplot(3,1,2)
plot(x(:,1), x(:,6), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_e\rm(m/s)')
grid on

subplot(3,1,3)
plot(x(:,1), x(:,7), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_d\rm(m/s)')
grid on

%-------------------------------------------------------------------------%
figure
subplot(3,1,1)
plot(x(:,1), x(:,8), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\psi(\circ)')
title('Attitude error')
grid on

subplot(3,1,2)
plot(x(:,1), x(:,9), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\theta(\circ)')
grid on

subplot(3,1,3)
plot(x(:,1), x(:,10), 'LineWidth',2)
set(gca, 'xlim', [x(1,1),x(end,1)])
xlabel('\itt\rm(s)')
ylabel('\delta\gamma(\circ)')
grid on

clearvars x