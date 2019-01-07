function plot_assembly_esti(assembly_esti, assembly_error, P)

t = assembly_esti(:,1);

P = P/pi*180;
P = P*3;

figure
subplot(3,1,1)
plot(t, assembly_esti(:,2), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,3)+assembly_error(1), 'Color','r', 'LineStyle','--');
plot(t, -P(:,3)+assembly_error(1), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\xi_\psi(\circ)')
title('安装误差角估计曲线')
grid on

subplot(3,1,2)
plot(t, assembly_esti(:,3), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,2)+assembly_error(2), 'Color','r', 'LineStyle','--');
plot(t, -P(:,2)+assembly_error(2), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\xi_\theta(\circ)')
grid on

subplot(3,1,3)
plot(t, assembly_esti(:,4), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,1)+assembly_error(3), 'Color','r', 'LineStyle','--');
plot(t, -P(:,1)+assembly_error(3), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\xi_\gamma(\circ)')
grid on

end