%% 1.Load data
imu = imu_delta;
traj = traj_m;

fs = 1/(imu(2,1)-imu(1,1));
[imu(:,2:4), ~] = imu_error(imu(:,2:4), fs, [0,0,0], 0.2/3600/180*pi, 0.04/3600/180*pi, 0.0015/3600/180*pi, 'delta');
[imu(:,5:7), ~] = imu_error(imu(:,5:7), fs, [0,0,0], 0.6*0.0098, 0.2*0.0098, 0.0055*0.0098, 'delta');

%% 2.Solve
earth_constant;
dt = (imu(2,1)-imu(1,1))*2;

n = size(imu,1)/2; %the number of inertial solving
nav = zeros(n,10); %[t, lat, lon, alt, vn, ve, vd, yaw, pitch, roll]

p = traj(1,2:4); %deg
v = traj(1,5:7); %m/s
att = traj(1,8:10); %deg
pva0 = [p, v, att]; %record initial value

att = att/180*pi; %rad
q = angle2quat(att(1), att(2), att(3));
Cnb = angle2dcm(att(1), att(2), att(3));
p(1:2) = p(1:2)/180*pi; %rad

p_1 = p;
p_2 = p;
v_1 = v;
v_2 = v;

for k=1:n
    kj = 2*k;
    dtheta1 = imu(kj-1,2:4);
    dtheta2 = imu(kj,  2:4);
    dv1 = imu(kj-1,5:7);
    dv2 = imu(kj,  5:7);
    
%     lat = (3*p_1(1)-p_2(1))/2; %extrapolation
%     h   = (3*p_1(3)-p_2(3))/2;
%     vn  = (3*v_1(1)-v_2(1))/2;
%     ve  = (3*v_1(2)-v_2(2))/2;
    lat = p(1); %last value
    h = p(3);
    vn = v(1);
    ve = v(2);
    
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    wien = [w*cos(lat), 0, -w*sin(lat)];
    wenn = [ve/Rn, -vn/Rm, -ve/Rn*tan(lat)];
    winn = wien+wenn;
    winb = (Cnb*winn')';
    
    %--velocity update--% v = traj(k+1,5:7);
    dtheta1_c = dtheta1 - winb*dt/2;
    dtheta2_c = dtheta2 - winb*dt/2;
    dvc = 1/2*cross(dtheta1_c,dv1)+7/6*cross(dtheta1_c,dv2)-1/6*cross(dtheta2_c,dv1)+1/2*cross(dtheta2_c,dv2);
    dv = (dv1+dv2+dvc)*Cnb - cross(2*wien+wenn,v)*dt + [0,0,9.8]*dt;
    if norm(dv)>1e-10
        v = v + dv;
    end
    
    %--position update--% p = [traj(k+1,2)/180*pi,traj(k+1,3)/180*pi,traj(k+1,4)]; 
    dp = (v+v_1)/2*dt;
    if norm(dp)>1e-10
        p(1) = p(1) + dp(1)/Rm;
        p(2) = p(2) + dp(2)/Rn*sec(lat);
        p(3) = p(3) - dp(3);
    end
%     p(3) = traj(k+1,4); %height locking
    
    %--attitude update--%
    Phi = dtheta1+dtheta2 + 2/3*cross(dtheta1,dtheta2);
    if norm(Phi)>1e-10
        phi = norm(Phi);
        qc = [cos(phi/2), Phi/phi*sin(phi/2)];
        q = quatmultiply(q, qc);
    end
    Cn0b = quat2dcm(q);
    Ce0n0 = dcmecef2ned(p_1(1)/pi*180,p_1(2)/pi*180);
    Ce0e = [cos(w*dt),sin(w*dt),0; -sin(w*dt),cos(w*dt),0; 0,0,1];
    Cen = dcmecef2ned(p(1)/pi*180,p(2)/pi*180);
    Cnb = Cn0b*Ce0n0*Ce0e'*Cen';
    q = dcm2quat(Cnb);
    
    v_2 = v_1;
    v_1 = v;
    p_2 = p_1;
    p_1 = p;
    
    %--store--%
    nav(k,1) = k*dt;
    nav(k,2:3) = p(1:2) /pi*180; %deg
    nav(k,4) = p(3);
    nav(k,5:7) = v;
    [r1,r2,r3] = dcm2angle(Cnb);
    nav(k,8:10) = [r1,r2,r3] /pi*180; %deg
end
nav = [[0,pva0]; nav];

%% 3.Plot
nav_error = plot_nav_error(traj, nav);