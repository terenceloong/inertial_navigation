sigma_gyro = 6/3600 /180*pi; %rad/s
sigma_acc = 0.06; %m/s^2
sigma_att = 0.01; %deg
sigma_v = 0.01; %m/s
sigma_lat = 0.1 /6378137; %rad
sigma_lon = 0.1 /6378137; %rad
sigma_h = 0.1; %m
gyro_bias = [0.3, -0.2, 0.1]*1; %deg/s
acc_bias = [0.02, 0.03, -0.04]*1; %m/s^2
assembly_error = [2.5, 2, -1]*1; %deg

%-------------------------------------------------------------------------%

Cf = angle2dcm(assembly_error(1)/180*pi, assembly_error(2)/180*pi, assembly_error(3)/180*pi);

n = size(imu,1);
for k=1:n
    imu(k,2:4) = (Cf*imu(k,2:4)')' + gyro_bias/180*pi; %rad/s
    imu(k,5:7) = (Cf*imu(k,5:7)')' + acc_bias; %m/s^2
end
imu(:,2:4) = imu(:,2:4) + randn(n,3)*sigma_gyro;
imu(:,5:7) = imu(:,5:7) + randn(n,3)*sigma_acc;

n = size(measure,1);
measure(:,2) = measure(:,2) + randn(n,1)*sigma_lat/pi*180; %deg
measure(:,3) = measure(:,3) + randn(n,1)*sigma_lon/pi*180; %deg
measure(:,4) = measure(:,4) + randn(n,1)*sigma_h;
measure(:,5:7) = measure(:,5:7) + randn(n,3)*sigma_v;
measure(:,8:10) = measure(:,8:10) + randn(n,3)*sigma_att; %deg

n = size(traj,1);
for k=1:n
    [r1,r2,r3] = dcm2angle(Cf*angle2dcm(traj(k,8)/180*pi,traj(k,9)/180*pi,traj(k,10)/180*pi));
    traj(k,8:10) = [r1,r2,r3] /pi*180; %deg
end