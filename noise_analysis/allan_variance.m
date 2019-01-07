function allan_variance(x, fs, option)
% Reference Matlab help, Inertial Sensor Noise Analysis Using Allan Variance
% N = 1, K = 0.1, --> B = 0.5

t0 = 1/fs;

if strcmp(option, 'omega')
    theta = cumsum(x,1)*t0;
else
    theta = cumsum(x,1);
end

m_max = 100;
L = length(theta);
m = round(10.^linspace(0,log10(floor((L-1)/2)),m_max))';
m = unique(m);
tau = m*t0;

sigma = zeros(length(tau),1);

for i=1:length(tau)
    mi = m(i);
    sigma(i) = sum((theta(1+2*mi:L) - 2*theta(1+mi:L-mi) + theta(1:L-2*mi)).^2);
end
sigma = sigma ./ (2*tau.^2.*(L-2*m));
sigma = sqrt(sigma);

% figure
% loglog(tau, sigma, 'LineWidth',2);
% title('Allan Deviation');
% xlabel('\tau');
% ylabel('\sigma(\tau)');
% grid on
% axis equal

%% Find N
slope = -0.5;
logtau = log10(tau);
logsigma = log10(sigma);
dlogsigma = diff(logsigma) ./ diff(logtau);
[~, i] = min(abs(dlogsigma - slope));
b = logsigma(i) - slope*logtau(i);
logN = slope*log10(1) + b;
N = 10^logN;
eval('N');

figure
loglog(tau, sigma, 'LineWidth',2);
hold on
loglog(tau, N./sqrt(tau), 'LineStyle','--', 'LineWidth',1);
loglog(1, N, 'bo');
text(1, N, '  N');
legend('\sigma', '\sigma_N');
title('Allan Deviation');
xlabel('\tau');
ylabel('\sigma(\tau)');
grid on
axis equal

%% Find B
slope = 0;
[~, i] = min(abs(dlogsigma - slope));
b = logsigma(i) - slope*logtau(i);
scfB = sqrt(2*log(2)/pi);
logB = b - log10(scfB);
B = 10^logB;
eval('B');
tauB = tau(i);

figure
loglog(tau, sigma, 'LineWidth',2);
hold on
loglog(tau, scfB*B*ones(length(tau),1), 'LineStyle','--', 'LineWidth',1);
loglog(tauB, scfB*B, 'bo');
text(tauB, scfB*B, '  0.664B');
legend('\sigma', '\sigma_B');
title('Allan Deviation');
xlabel('\tau');
ylabel('\sigma(\tau)');
grid on
axis equal

%% Find K
slope = 0.5;
[~, i] = min(abs(dlogsigma - slope));
b = logsigma(i) - slope*logtau(i);
logK = slope*log10(3) + b;
K = 10^logK;
eval('K');

figure
loglog(tau, sigma, 'LineWidth',2);
hold on
loglog(tau, K.*sqrt(tau/3), 'LineStyle','--', 'LineWidth',1);
loglog(3, K, 'bo');
text(3, K, '  K');
legend('\sigma', '\sigma_K');
title('Allan Deviation');
xlabel('\tau');
ylabel('\sigma(\tau)');
grid on
axis equal

%% Plot
figure
loglog(tau, sigma, 'LineWidth',2);
hold on
loglog(tau, N./sqrt(tau), 'LineStyle','--', 'LineWidth',1);
loglog(tau, scfB*B*ones(length(tau),1), 'LineStyle','--', 'LineWidth',1);
loglog(tau, K.*sqrt(tau/3), 'LineStyle','--', 'LineWidth',1);
loglog(1, N, 'bo');
text(1, N, '  N');
loglog(tauB, scfB*B, 'bo');
text(tauB, scfB*B, '  0.664B');
loglog(3, K, 'bo');
text(3, K, '  K');
legend('\sigma', '\sigma_N', '\sigma_B', '\sigma_K');
title('Allan Deviation');
xlabel('\tau');
ylabel('\sigma(\tau)');
grid on
axis equal

end