function qvp = ins_qvp_n3(qvp, dt, imu0, imu1, imu2)

global a f w
q = qvp(1:4);
v = qvp(5:7);
v0 = v;
lat = qvp(8);
lon = qvp(9);
h = qvp(10);
Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
Cnb = quat2dcm(q');
wien = [w*cos(lat); 0; -w*sin(lat)];
wenn = [v(2)/Rn; -v(1)/Rm; -v(2)/Rn*tan(lat)];
winb = Cnb*(wien+wenn);  
wnbb0 = imu0(1:3) - winb;
wnbb1 = imu1(1:3) - winb;
wnbb2 = imu2(1:3) - winb;
acc0 = imu0(4:6);
acc1 = imu1(4:6);
acc2 = imu2(4:6);
dtheta1 = (wnbb0+wnbb1)*dt/4;
dtheta2 = (wnbb1+wnbb2)*dt/4;
dv1 = (acc0+acc1)*dt/4;
dv2 = (acc1+acc2)*dt/4;
q = RK4(@fun_dq, q, dt, wnbb0,wnbb1,wnbb2);
q = quatnormalize(q')';
dvc = 0.5*cross(dtheta1,dv1) + 7/6*cross(dtheta1,dv2) - 1/6*cross(dtheta2,dv1) + 0.5*cross(dtheta2,dv2);
dv = Cnb'*(dv1+dv2+dvc) - dt*cross((2*wien+wenn),v) + dt*[0;0;9.8];
% if norm(dv)>1e-10
    v = v + dv;
% end
dp = (v0+v)/2*dt;
% if norm(dp)>1e-10
    lon = lon + dp(2)/(Rn*cos(lat));
    lat = lat + dp(1)/Rm;
    h = h - dp(3);
% end
qvp = [q; v ; lat; lon; h];
    
end

function dq = fun_dq(q, w)
dq = 0.5*[ 0,   -w(1), -w(2), -w(3);
          w(1),   0,    w(3), -w(2);
          w(2), -w(3),   0,    w(1);
          w(3),  w(2), -w(1),   0 ]*q;
end