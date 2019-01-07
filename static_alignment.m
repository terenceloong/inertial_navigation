% 零偏不要反馈，因为零偏反馈时需要准确的姿态角，对于高精度的情况，姿态角总会有误差，反馈会带来稳态误差
% 零偏状态量选取为导航系下零偏，东向陀螺零偏不可估，所以不加，如果加了，航向估计不动，航向误差全归到东向陀螺零偏上去，怎么分目前不知道
% 北向陀螺零偏估计快，天向陀螺零偏估计慢
% 天向陀螺零偏要等航向角稳定后才能出来
% 陀螺零偏大了，对应的P也加大

% clear; clc;

%% 1.参数设置
global lat w g
lat = 45; %deg 纬度
w = 15; %deg/h 地球自转角速度
g = 9.8; %m/s^2 重力加速度
att = [0, 0, 0]; %deg 真实姿态（航向，俯仰，滚转）
datt = [-0.6, 0.1, -0.1]*1; %deg 惯导解算初始姿态误差（航向，俯仰，滚转）
bias_gyro = [1; 0.1; -1]; %deg/h 陀螺仪零偏
bias_acc = [1; 1; 0.5]*1e-3; %g 加速度计零偏
dt = 0.01; %s 惯导采样周期

%% 2.生成惯导数据
Cnb = angle2dcm(att(1)/180*pi, att(2)/180*pi ,att(3)/180*pi);
wn = [cosd(lat);0;-sind(lat)]*w; %deg/h
fn = [0;0;-1]; %g
wb = Cnb*wn + bias_gyro; %deg/h
fb = Cnb*fn + bias_acc; %g

sigma_gyro = 0.2; %deg/h 陀螺仪角速度噪声标准差
sigma_acc = 0.1e-3; %g 加速度计噪声标准差
sigma_v = 1e-3; %m/s 速度噪声标准差，滤波器用

t = 0:dt:1200;
n = length(t);
imu = [ones(n,1)*wb'+randn(n,3)*sigma_gyro, ones(n,1)*fb'+randn(n,3)*sigma_acc]; %deg/h, g

%% 3.惯导解算初值计算
n = (n-1)/2;
dt = dt*2;

nav = zeros(n,6); %存储导航结果

q = angle2quat((att(1)+datt(1))/180*pi, (att(2)+datt(2))/180*pi, (att(3)+datt(3))/180*pi);
av = [q'; [0;0;0]];

%% 4.滤波器参数
system_model = 'model_4';
switch system_model
    case 'model_1'
        N = 6; %滤波器维数
    case 'model_2'
        N = 8;
    case 'model_3'
        N = 9;
    case 'model_4'
        N = 9;
    case 'model_5'
        N = 7;
end
X = zeros(N,1);
H = zeros(3,N);
H(1:3,4:6) = eye(3);
A = zeros(N);
A(1:3,1:3) = -antisym([cosd(lat),0,-sind(lat)]) *w/3600/180*pi;
A(4:6,1:3) = antisym([0,0,-g]);
switch system_model
    case 'model_1'
        P = diag([[1,1,1]*0.5/180*pi, [1,1,1]*0.1].^2);
        Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1, [1,1,1]*(sigma_acc*g)*0.707*1].^2) *dt^2;
	case 'model_2'
        A(1:3,7:8) = -[1,0; 0,0; 0,1];
        P = diag([[1,1,1]*0.5/180*pi, [1,1,1]*0.1, [1,1]*1/3600/180*pi].^2);
        Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1, [1,1,1]*(sigma_acc*g)*0.707*1, [0,0]*(1/3600/180*pi)].^2) *dt^2;
    case 'model_3'
        A(1:3,7:9) = -eye(3);
        P = diag([[1,1,1]*0.5/180*pi, [1,1,1]*0.1, [1,1,1]*1/3600/180*pi].^2);
        Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1, [1,1,1]*(sigma_acc*g)*0.707*1, [0,0,0]*(1/3600/180*pi)].^2) *dt^2;
    case 'model_4'
        A(1:3,7:8) = -[1,0; 0,0; 0,1];
        A(6,9) = 1;
        P = diag([[1,1,1]*0.5/180*pi, [1,1,1]*0.1, [1,1]*1/3600/180*pi, 1*1e-2].^2);
        Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1, [1,1,1]*(sigma_acc*g)*0.707*1, [0,0]*(1/3600/180*pi), 0].^2) *dt^2;
    case 'model_5'
        A(1:3,7) = -[1; 0; 0];
        P = diag([[1,1,1]*0.5/180*pi, [1,1,1]*0.1, 1*1/3600/180*pi].^2);
        Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1, [1,1,1]*(sigma_acc*g)*0.707*1, 0*(1/3600/180*pi)].^2) *dt^2;
end
Phi = eye(N)+A*dt+(A*dt)^2/2;
P0 = sqrt(diag(P))';
R = diag([1,1,1]*sigma_v^2);

bias_esti = zeros(n,6); %存零偏估计
filter_P = zeros(n,N); %存滤波器P阵

%% 5.算法执行
for k=1:n
    kj = 2*k+1;
    gyro0 = imu(kj-2, 1:3)' /3600/180*pi; %rad/s
    gyro1 = imu(kj-1, 1:3)' /3600/180*pi;
    gyro2 = imu(kj  , 1:3)' /3600/180*pi;
    acc0  = imu(kj-2, 4:6)' *g; %m/s^2
    acc1  = imu(kj-1, 4:6)' *g;
    acc2  = imu(kj  , 4:6)' *g;
    
    av = RK4(@fun_dx, av, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
    
    %---------------------------------------------------------------------%
    Z = av(5:7);% + randn(3,1)*sigma_v;
    [X, P] = kalman_filter(Phi, X, P, Q, H, Z, R);
    filter_P(k,:) = sqrt(diag(P))';
    
    %---------adjust----------%
    if norm(X(1:3))>0
        phi = norm(X(1:3));
        qc = [cos(phi/2), X(1:3)'/phi*sin(phi/2)];
        av(1:4) = quatmultiply(qc, av(1:4)')';
    end
    av(5:7) = av(5:7) - X(4:6);
    X(1:6) = zeros(6,1);
    switch system_model
        case 'model_2'
            bias_esti(k,1) = X(7) /pi*180*3600; %deg/h
            bias_esti(k,3) = X(8) /pi*180*3600; %deg/h
    	case 'model_3'
            bias_esti(k,1) = X(7) /pi*180*3600; %deg/h
            bias_esti(k,2) = X(8) /pi*180*3600; %deg/h
            bias_esti(k,3) = X(9) /pi*180*3600; %deg/h
        case 'model_4'
            bias_esti(k,1) = X(7) /pi*180*3600; %deg/h
            bias_esti(k,3) = X(8) /pi*180*3600; %deg/h
            bias_esti(k,6) = X(9) /g*1000; %mg
        case 'model_5'
            bias_esti(k,1) = X(7) /pi*180*3600; %deg/h
    end
    %---------------------------------------------------------------------%
    
    [r1,r2,r3] = quat2angle(av(1:4)');
    nav(k,1:3) = [r1,r2,r3] /pi*180; %deg
    nav(k,4:6) = av(5:7)';
end
nav = [[(att+datt),0,0,0];nav];
filter_P = [P0; filter_P];
bias_esti = [zeros(1,6); bias_esti];

%% 6.画图
figure
subplot(3,2,1)
plot(t(1:2:end), nav(:,1)-att(1))
hold on
axis manual
plot(t(1:2:end),  filter_P(:,3)/pi*180*3, 'Color','r', 'LineStyle','--');
plot(t(1:2:end), -filter_P(:,3)/pi*180*3, 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\psi(\circ)')
grid on

subplot(3,2,3)
plot(t(1:2:end), nav(:,2)-att(2))
hold on
axis manual
plot(t(1:2:end),  filter_P(:,2)/pi*180*3, 'Color','r', 'LineStyle','--');
plot(t(1:2:end), -filter_P(:,2)/pi*180*3, 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\theta(\circ)')
grid on

subplot(3,2,5)
plot(t(1:2:end), nav(:,3)-att(3))
hold on
axis manual
plot(t(1:2:end),  filter_P(:,1)/pi*180*3, 'Color','r', 'LineStyle','--');
plot(t(1:2:end), -filter_P(:,1)/pi*180*3, 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\gamma(\circ)')
grid on

subplot(3,2,2)
plot(t(1:2:end), nav(:,4))
hold on
axis manual
plot(t(1:2:end),  filter_P(:,4)*3, 'Color','r', 'LineStyle','--');
plot(t(1:2:end), -filter_P(:,4)*3, 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\itv_n\rm(m/s)')
grid on

subplot(3,2,4)
plot(t(1:2:end), nav(:,5))
hold on
axis manual
plot(t(1:2:end),  filter_P(:,5)*3, 'Color','r', 'LineStyle','--');
plot(t(1:2:end), -filter_P(:,5)*3, 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\itv_e\rm(m/s)')
grid on

subplot(3,2,6)
plot(t(1:2:end), nav(:,6))
hold on
axis manual
plot(t(1:2:end),  filter_P(:,6)*3, 'Color','r', 'LineStyle','--');
plot(t(1:2:end), -filter_P(:,6)*3, 'Color','r', 'LineStyle','--');
set(gca, 'xlim', [t(1),t(end)])
ylabel('\delta\itv_d\rm(m/s)')
grid on

if strcmp(system_model,'model_2') || strcmp(system_model,'model_4')
    figure
    subplot(2,1,1)
    plot(t(1:2:end), bias_esti(:,1))
    hold on
    axis manual
    plot(t(1:2:end),  filter_P(:,7)/pi*180*3600*3+bias_gyro(1), 'Color','r', 'LineStyle','--');
    plot(t(1:2:end), -filter_P(:,7)/pi*180*3600*3+bias_gyro(1), 'Color','r', 'LineStyle','--');
    set(gca, 'xlim', [t(1),t(end)])
    xlabel('\itt\rm(s)')
    ylabel('\it\epsilon_x\rm(\circ/h)')
    grid on
    
    subplot(2,1,2)
    plot(t(1:2:end), bias_esti(:,3))
    hold on
    axis manual
    plot(t(1:2:end),  filter_P(:,8)/pi*180*3600*3+bias_gyro(3), 'Color','r', 'LineStyle','--');
    plot(t(1:2:end), -filter_P(:,8)/pi*180*3600*3+bias_gyro(3), 'Color','r', 'LineStyle','--');
    set(gca, 'xlim', [t(1),t(end)])
    xlabel('\itt\rm(s)')
    ylabel('\it\epsilon_z\rm(\circ/h)')
    grid on
end

if strcmp(system_model,'model_4')
    figure
    plot(t(1:2:end), bias_esti(:,6))
    hold on
    axis manual
    plot(t(1:2:end),  filter_P(:,9)/g*1000*3+bias_acc(3)*1000, 'Color','r', 'LineStyle','--');
    plot(t(1:2:end), -filter_P(:,9)/g*1000*3+bias_acc(3)*1000, 'Color','r', 'LineStyle','--');
    set(gca, 'xlim', [t(1),t(end)])
    xlabel('\itt\rm(s)')
    ylabel('\it\nabla_z\rm(mg)')
    grid on
end

%%
% figure
% set(gcf,'position',[680,216,560,760])
% subplot(3,1,1)
% plot(t(1:2:end), nav(:,1)-att(1), 'LineWidth',2)
% hold on
% axis manual
% plot(t(1:2:end),  filter_P(:,3)/pi*180*3, 'Color','r', 'LineStyle','--');
% plot(t(1:2:end), -filter_P(:,3)/pi*180*3, 'Color','r', 'LineStyle','--');
% set(gca, 'xlim', [t(1),t(end)])
% xlabel('\itt\rm(s)')
% ylabel('\delta\psi(\circ)')
% title('航向角估计误差曲线')
% grid on
% 
% subplot(3,1,2)
% plot(t(1:2:end), nav(:,2)-att(2), 'LineWidth',2)
% hold on
% axis manual
% plot(t(1:2:end),  filter_P(:,2)/pi*180*3, 'Color','r', 'LineStyle','--');
% plot(t(1:2:end), -filter_P(:,2)/pi*180*3, 'Color','r', 'LineStyle','--');
% set(gca, 'xlim', [t(1),t(end)])
% xlabel('\itt\rm(s)')
% ylabel('\delta\theta(\circ)')
% title('俯仰角估计误差曲线')
% grid on
% 
% subplot(3,1,3)
% plot(t(1:2:end), nav(:,3)-att(3), 'LineWidth',2)
% hold on
% axis manual
% plot(t(1:2:end),  filter_P(:,1)/pi*180*3, 'Color','r', 'LineStyle','--');
% plot(t(1:2:end), -filter_P(:,1)/pi*180*3, 'Color','r', 'LineStyle','--');
% set(gca, 'xlim', [t(1),t(end)])
% xlabel('\itt\rm(s)')
% ylabel('\delta\gamma(\circ)')
% title('滚转角估计误差曲线')
% grid on
% 
% figure
% set(gcf,'position',[680,216,560,760])
% subplot(3,1,1)
% plot(t(1:2:end), bias_esti(:,1), 'LineWidth',2)
% hold on
% axis manual
% plot(t(1:2:end),  filter_P(:,7)/pi*180*3600*3+bias_gyro(1), 'Color','r', 'LineStyle','--');
% plot(t(1:2:end), -filter_P(:,7)/pi*180*3600*3+bias_gyro(1), 'Color','r', 'LineStyle','--');
% set(gca, 'xlim', [t(1),t(end)])
% xlabel('\itt\rm(s)')
% ylabel('\it\epsilon_x\rm(\circ/h)')
% title('北向陀螺仪零偏估计曲线')
% grid on
% 
% subplot(3,1,2)
% plot(t(1:2:end), bias_esti(:,3), 'LineWidth',2)
% hold on
% axis manual
% plot(t(1:2:end),  filter_P(:,8)/pi*180*3600*3+bias_gyro(3), 'Color','r', 'LineStyle','--');
% plot(t(1:2:end), -filter_P(:,8)/pi*180*3600*3+bias_gyro(3), 'Color','r', 'LineStyle','--');
% set(gca, 'xlim', [t(1),t(end)])
% xlabel('\itt\rm(s)')
% ylabel('\it\epsilon_z\rm(\circ/h)')
% title('地向陀螺仪零偏估计曲线')
% grid on
% 
% subplot(3,1,3)
% plot(t(1:2:end), bias_esti(:,6), 'LineWidth',2)
% hold on
% axis manual
% plot(t(1:2:end),  filter_P(:,9)/g*1000*3+bias_acc(3)*1000, 'Color','r', 'LineStyle','--');
% plot(t(1:2:end), -filter_P(:,9)/g*1000*3+bias_acc(3)*1000, 'Color','r', 'LineStyle','--');
% set(gca, 'xlim', [t(1),t(end)])
% xlabel('\itt\rm(s)')
% ylabel('\it\nabla_z\rm(mg)')
% title('地向加速度计零偏估计曲线')
% grid on

%% Function
function dx = fun_dx(x, u)
% x = [q0;q1;q2;q3; vx;vy;vz]
% u = [wibbx;wibby;wibbz; fbx;fby;fbz]
    global lat w g
    q = x(1:4);
    wibb = u(1:3);
    fb = u(4:6);

    Cnb = quat2dcm(q');
    winb = Cnb*[cosd(lat);0;-sind(lat)]*w/180*pi/3600;
    wnbb = wibb - winb;

    dq = 0.5*[ 0,   -wnbb(1), -wnbb(2), -wnbb(3);
              wnbb(1),   0,    wnbb(3), -wnbb(2);
              wnbb(2), -wnbb(3),   0,    wnbb(1);
              wnbb(3),  wnbb(2), -wnbb(1),   0 ]*q;
    fn = Cnb'*fb + [0;0;g];
    dx =[dq; fn];
end