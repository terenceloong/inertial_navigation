clear;clc;

T = 80;
traj0 = traj_design(0);
p0 = [30, 120, 300]; %deg, [lat,lon,h]

%---------------------------------------------------------------------%
v0 = traj0(4:6);
att0 = traj0(1:3);
Cnb = angle2dcm(att0(1),att0(2),att0(3));
vb0 = v0*Cnb';
att0 = [att0(3),att0(2),att0(1)];

a = 6378137;
f = 1/298.257223563;
w = 7.292115e-5;
lat = p0(1)/180*pi;
h = p0(3);
v = v0;
Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
wien = [w*cos(lat), 0, -w*sin(lat)];
wenn = [v(2)/Rn, -v(1)/Rm, -v(2)/Rn*tan(lat)];
f_harmful = cross(2*wien+wenn,v);

clearvars -except T p0 vb0 att0 f_harmful
%---------------------------------------------------------------------%

dt_scope = 0.1;
dt_solve = 0.01;
dt_imu_rate = 0.01;
dt_traj = 0.02;