function error = plot_nav_error(traj, nav, P)

t = traj(:,1);

error = nav(:,2:end) - traj(:,2:end);
error(:,1:2) = error(:,1:2)/180*pi*6378137; %lat and lon, deg to m
for k=1:size(error,1)
    if error(k,7)>300
        error(k,7) = error(k,7)-360;
    elseif error(k,7)<-300
        error(k,7) = error(k,7)+360;
    end
    if error(k,9)>300
        error(k,9) = error(k,9)-360;
    elseif error(k,9)<-300
        error(k,9) = error(k,9)+360;
    end
end

if nargin==3
    P(:,1:3) = P(:,1:3)/pi*180;
    P(:,7:8) = P(:,7:8)*6378137;
    P = P*3; %3 sigma
end

figure
% subplot(3,3,1)
subplot(3,1,1)
plot(t, error(:,1), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,7), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,7), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itL\rm(m)')
title('Î»ÖÃÎó²îÇúÏß')
grid on

% subplot(3,3,4)
subplot(3,1,2)
plot(t, error(:,2), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,8), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,8), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\lambda(m)')
grid on

% subplot(3,3,7)
subplot(3,1,3)
plot(t, error(:,3), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,9), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,9), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\ith\rm(m)')
grid on

% subplot(3,3,2)
figure
subplot(3,1,1)
plot(t, error(:,4), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,4), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,4), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_n\rm(m/s)')
title('ËÙ¶ÈÎó²îÇúÏß')
grid on

% subplot(3,3,5)
subplot(3,1,2)
plot(t, error(:,5), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,5), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,5), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_e\rm(m/s)')
grid on

% subplot(3,3,8)
subplot(3,1,3)
plot(t, error(:,6), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,6), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,6), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\itv_d\rm(m/s)')
grid on

% subplot(3,3,3)
figure
subplot(3,1,1)
plot(t, error(:,7), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,3), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,3), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\psi(\circ)')
title('×ËÌ¬Îó²îÇúÏß')
grid on

% subplot(3,3,6)
subplot(3,1,2)
plot(t, error(:,8), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,2), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,2), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\theta(\circ)')
grid on

% subplot(3,3,9)
subplot(3,1,3)
plot(t, error(:,9), 'LineWidth',2)
if nargin==3
    hold on
    axis manual
    plot(t,  P(:,1), 'Color','r', 'LineStyle','--');
    plot(t, -P(:,1), 'Color','r', 'LineStyle','--');
end
set(gca, 'xlim', [t(1),t(end)])
xlabel('\itt\rm(s)')
ylabel('\delta\gamma(\circ)')
grid on

error = [traj(:,1),error];

end