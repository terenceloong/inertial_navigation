function plot_gyro_esti(bias_esti, gyro_bias, P)

t = bias_esti(:,1);

P = P/pi*180;
P = P*3;

figure
subplot(3,1,1)
plot(t, bias_esti(:,2), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,1)+gyro_bias(1), 'Color','r', 'LineStyle','--');
plot(t, -P(:,1)+gyro_bias(1), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\epsilon_x\rm(\circ/s)')
title('ÍÓÂÝÒÇÁãÆ«¹À¼ÆÇúÏß')
grid on

subplot(3,1,2)
plot(t, bias_esti(:,3), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,2)+gyro_bias(2), 'Color','r', 'LineStyle','--');
plot(t, -P(:,2)+gyro_bias(2), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\epsilon_y\rm(\circ/s)')
grid on

subplot(3,1,3)
plot(t, bias_esti(:,4), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,3)+gyro_bias(3), 'Color','r', 'LineStyle','--');
plot(t, -P(:,3)+gyro_bias(3), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\epsilon_z\rm(\circ/s)')
grid on

end