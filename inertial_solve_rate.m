%% 1.Load data
imu = imu_rate(:,1:7);
traj = traj_m;

% fs = 1/(imu(2,1)-imu(1,1));
% [imu(:,2:4), ~] = imu_error(imu(:,2:4), fs, [0,0,0], 0.2/3600/180*pi, 0.04/3600/180*pi, 0.0015/3600/180*pi, 'rate');
% [imu(:,5:7), ~] = imu_error(imu(:,5:7), fs, [0,0,0], 0.6*0.0098, 0.2*0.0098, 0.0055*0.0098, 'rate');

%% 2.Solve
earth_constant;
dt = (imu(2,1)-imu(1,1))*2;

n = (size(imu,1)-1)/2; %the number of inertial solving
nav = zeros(n,10); %[t, lat, lon, alt, vn, ve, vd, yaw, pitch, roll]

p = traj(1,2:4); %deg
v = traj(1,5:7); %m/s
att = traj(1,8:10); %deg
pva0 = [p, v, att]; %record initial value

p(1:2) = p(1:2)/180*pi; %rad
att = att/180*pi; %rad
q = angle2quat(att(1), att(2), att(3));
qvp = [q, v, p]'; %column vector

for k=1:n
    kj = 2*k+1;
    gyro0 = imu(kj-2, 2:4)'; %rad/s
    gyro1 = imu(kj-1, 2:4)';
    gyro2 = imu(kj  , 2:4)';
    acc0  = imu(kj-2, 5:7)'; %m/s^2
    acc1  = imu(kj-1, 5:7)';
    acc2  = imu(kj  , 5:7)';
    
%     qvp = RK4(@ins_qvp_n, qvp, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
%     qvp(1:4) = quatnormalize(qvp(1:4)')'; %quaternion normalization
    qvp = ins_qvp_n3(qvp, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
    
%     qvp(10) = traj(k+1,4); %height locking
    
    nav(k,1) = k*dt;
    nav(k,2:3) = qvp(8:9)' /pi*180; %deg
    nav(k,4) = qvp(10); %m
    nav(k,5:7) = qvp(5:7)'; %m/s
    [r1,r2,r3] = quat2angle(qvp(1:4)');
    nav(k,8:10) = [r1,r2,r3] /pi*180; %deg
end
nav = [[0,pva0]; nav];

%% 3.Plot
plot_nav_error(traj, nav);