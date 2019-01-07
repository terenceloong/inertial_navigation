function plot_acc_esti(bias_esti, acc_bias, P)

t = bias_esti(:,1);

bias_esti(:,5:7) = bias_esti(:,5:7)*100;
acc_bias = acc_bias*100;

P = P*100;
P = P*3;

figure
subplot(3,1,1)
plot(t, bias_esti(:,5), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,1)+acc_bias(1), 'Color','r', 'LineStyle','--');
plot(t, -P(:,1)+acc_bias(1), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\nabla_x\rm(mg)')
title('加速度计零偏估计曲线')
grid on

subplot(3,1,2)
plot(t, bias_esti(:,6), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,2)+acc_bias(2), 'Color','r', 'LineStyle','--');
plot(t, -P(:,2)+acc_bias(2), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\nabla_y\rm(mg)')
grid on

subplot(3,1,3)
plot(t, bias_esti(:,7), 'LineWidth',2)
hold on
axis manual
plot(t,  P(:,3)+acc_bias(3), 'Color','r', 'LineStyle','--');
plot(t, -P(:,3)+acc_bias(3), 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\it\nabla_z\rm(mg)')
grid on

end