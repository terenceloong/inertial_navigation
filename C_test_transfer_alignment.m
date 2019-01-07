%% 1.1 Precise data, add noise
imu = imu_rate(:,1:7); %[t, wx,wy,wz, fx,fy,fz]
measure = traj_m; %[t, lat,lon,h, vn,ve,vd, psi,theta,gamma], deg
traj = traj_m;

add_noise;

%% 1.2 Existing data
% imu
% measure
% traj

% sigma_gyro = 0.1/3600 /180*pi; %rad/s
% sigma_acc = 0.001; %m/s^2
% sigma_att = 0.01; %deg
% sigma_v = 0.01; %m/s
% sigma_lat = 0.1 /6378137; %rad
% sigma_lon = 0.1 /6378137; %rad
% sigma_h = 0.1; %m
% gyro_bias = [0, 0, 0]; %deg/s
% acc_bias = [0, 0, 0]; %m/s^2
% assembly_error = [0, 0, 0]; %deg

%% 2.Initial value
p0_error   = [0, 0, 0]; %deg
v0_error   = [0, 0, 0]; %m/s
att0_error = [2.5, 2, -1]*0; %deg

P0_phi = 1 /180*pi; %rad
P0_dv = 1; %m/s
P0_dlat = 5 /6378137; %rad
P0_dlon = 5 /6378137; %rad
P0_dh = 5; %m
P0_xi = 1 /180*pi; %rad
P0_e = 1 /180*pi; %rad/s
P0_d = 0.02; %m/s^2

dt = imu(3,1) - imu(1,1);

system_model = 'model_5';

% P,Q,R is vector, they are the diagonal element
switch system_model
    case 'model_1' %姿态、速度、位置
        P = [[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh]].^2;
        Q = [[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2].^2 * dt^2;
        R = [[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h]].^2;
    case 'model_3' %姿态、速度、位置、安装误差角、陀螺仪零偏
        P = [[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh], [1,1,1]*P0_xi, [1,1,1]*P0_e].^2;
        Q = [[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2, [1,1,1]*0, [1,1,1]*0].^2 * dt^2;
        R = [[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h], [1,1,1]*0.001].^2;
    case 'model_5' %姿态、速度、位置、安装误差角、陀螺仪零偏、加速度计零偏
        P = [[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh], [1,1,1]*P0_xi, [1,1,1]*P0_e, [1,1,1]*P0_d].^2;
        Q = [[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2, [1,1,1]*0, [1,1,1]*0, [1,1,1]*0].^2 * dt^2;
        R = [[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h], [1,1,1]*0.001].^2;
end

p = traj(1,2:4) + p0_error; %deg
v = traj(1,5:7) + v0_error; %m/s
att = traj(1,8:10) + att0_error; %deg
p(1:2) = p(1:2)/180*pi; %rad
att = att/180*pi; %rad
q = angle2quat(att(1), att(2), att(3));
qvp = [q, v, p]';

%% 3.Transfer alignment
[nav, bias_esti, assembly_esti, filter_P] = test_transfer_alignment(imu, measure, qvp, P, Q, R);

%% 4.Plot navigation error
plot_nav_error(traj, nav, filter_P(:,1:9));
if strcmp(system_model,'model_2')
    plot_assembly_esti(assembly_esti, assembly_error, filter_P(:,10:12));
end
if strcmp(system_model,'model_3')
    plot_assembly_esti(assembly_esti, assembly_error, filter_P(:,10:12));
    plot_gyro_esti(bias_esti, gyro_bias, filter_P(:,13:15));
end
if strcmp(system_model,'model_4')
    plot_gyro_esti(bias_esti, gyro_bias, filter_P(:,10:12));
end
if strcmp(system_model,'model_5')
    plot_assembly_esti(assembly_esti, assembly_error, filter_P(:,10:12));
    plot_gyro_esti(bias_esti, gyro_bias, filter_P(:,13:15));
    plot_acc_esti(bias_esti, acc_bias, filter_P(:,16:18));
end