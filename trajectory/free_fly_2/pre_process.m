%% 1.Load data
clear;clc;

path = './free_fly_2';
data = load('Data.txt');
load('Coef.mat');

%% 2.Calculate velocity
v = [[0,0,0]; diff(data(:,2:4))./diff(data(:,1))];

a = 6378137;
f = 1/298.257223563;
for k=1:size(v,1)
    lat = data(k,2)/180*pi; %rad
    h = data(k,4);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    v(k,1) = v(k,1)/180*pi*Rm;
    v(k,2) = v(k,2)/180*pi*Rn*cos(lat);
end

%% 3.Fliter
N = (length(Coef)-1)/2;

in = [zeros(N,3); v]; %add N of zero
vf = zeros(size(in,1)-2*N,3);
for k=N+1:size(in,1)-N
    vf(k-N,1:3) = sum(in(k-N:k+N,1:3).*Coef',1);
end

v(end-N+1:end,:) = [];

t = data(1:size(v,1),1)-data(1,1);

figure
plot(t, v(:,1))
hold on
plot(t, vf(:,1), 'LineWidth',2)

figure
plot(t, v(:,2))
hold on
plot(t, vf(:,2), 'LineWidth',2)

figure
plot(t, v(:,3))
hold on
plot(t, vf(:,3), 'LineWidth',2)

v = vf;

%% 4.Generate cmd
data = data(1:size(v,1),:);
data(:,1) = data(:,1)-data(1,1);
n = size(data,1);

cmd = zeros(n,8); %[t, vn, ve, vd, q0,q1,q2,q3]
cmd(:,1) = data(:,1); %t
cmd(:,2:4) = v; %velocity
cmd(:,4) = -cmd(:,4);
for k=1:n
    r1 = round(data(k,5),5)/180*pi;
    r2 = round(data(k,6),5)/180*pi;
    r3 = round(data(k,7),5)/180*pi;
    cmd(k,5:8) = angle2quat(r1, r2, r3); %q
end
for k=2:n %make quaternion continue
    if abs(cmd(k,5)-cmd(k-1,5))>1
        cmd(k:end,5:8) = -cmd(k:end,5:8);
    end
end
pos0 = data(1,2:4); %deg
att0 = round(data(1,[7,6,5]),5) /180*pi; %rad, [roll, pitch, yaw]

%% 5.Set time parameter
dt_scope = 0.1;
dt_solve = 1e-3;
dt_imu_rate = 10e-3;
dt_imu_delta = 10e-3;
dt_traj = 20e-3;
Pv = 0.3;

%% 6.Save
save('cmd.mat', 'data','cmd','pos0','att0','dt_scope','dt_solve','dt_imu_rate','dt_imu_delta','dt_traj','path','Pv');
clear