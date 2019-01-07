% imu = imu_rate(:,1:7);
imu = imu_delta;
traj = traj_m;

%% 1.Initial value
p = traj(1,2:4);
p(1:2) = p(1:2)/180*pi;
v = traj(1,5:7);
att = traj(1,8:10)/180*pi;
q = angle2quat(att(1), att(2), att(3));
avp = [q, v, p]';

%% 2.Inertial solve
nav = test_inertial_solve(imu, avp, traj(:,4));

%% 3.Plot navigation error
plot_nav_error(traj, nav);