%% 1.Load data
imu = imu_rate(:,1:7); %[t, wx,wy,wz, fx,fy,fz]
measure = traj_m; %[t, lat,lon,h, vn,ve,vd, psi,theta,gamma], deg
traj = traj_m;

add_noise;

%% 2.Initial value
earth_constant;

p0_error   = [0, 0, 0]; %deg
v0_error   = [0, 0, 0]; %m/s
att0_error = [2.5, 2, -1]*0; %deg
system_model = 'model_5';

P0_phi = 1 /180*pi; %rad
P0_dv = 1; %m/s
P0_dlat = 5 /6378137; %rad
P0_dlon = 5 /6378137; %rad
P0_dh = 5; %m
P0_xi = 1 /180*pi; %rad
P0_e = 1 /180*pi; %rad/s
P0_d = 0.02; %m/s^2

dt = imu(3,1) - imu(1,1);
switch system_model
    case 'model_1' %×ËÌ¬¡¢ËÙ¶È¡¢Î»ÖÃ
        N = 9;
        H = zeros(6,N);
        H(1:3,4:6) = eye(3);
        H(4:6,7:9) = eye(3);
        P = diag([[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh]])^2;
        Q = diag([[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2])^2 * dt^2;
        R = diag([[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h]])^2;
    case 'model_2' %×ËÌ¬¡¢ËÙ¶È¡¢Î»ÖÃ¡¢°²×°Îó²î½Ç
        N = 12;
        H = zeros(9,N);
        H(1:3,4:6) = eye(3);
        H(4:6,7:9) = eye(3);
        H(7:9,1:3) = eye(3);
        P = diag([[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh], [1,1,1]*P0_xi])^2;
        Q = diag([[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2, [1,1,1]*0])^2 * dt^2;
        R = diag([[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h], [1,1,1]*0.001])^2;
    case 'model_3' %×ËÌ¬¡¢ËÙ¶È¡¢Î»ÖÃ¡¢°²×°Îó²î½Ç¡¢ÍÓÂÝÒÇÁãÆ«
        N = 15;
        H = zeros(9,N);
        H(1:3,4:6) = eye(3);
        H(4:6,7:9) = eye(3);
        H(7:9,1:3) = eye(3);
        P = diag([[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh], [1,1,1]*P0_xi, [1,1,1]*P0_e])^2;
        Q = diag([[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2, [1,1,1]*0, [1,1,1]*0])^2 * dt^2;
        R = diag([[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h], [1,1,1]*0.001])^2;
    case 'model_4' %×ËÌ¬¡¢ËÙ¶È¡¢Î»ÖÃ¡¢ÍÓÂÝÒÇÁãÆ«
        N = 12;
        H = zeros(6,N);
        H(1:3,4:6) = eye(3);
        H(4:6,7:9) = eye(3);
        P = diag([[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh], [1,1,1]*P0_e])^2;
        Q = diag([[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2, [1,1,1]*0])^2 * dt^2;
        R = diag([[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h]])^2;
    case 'model_5' %×ËÌ¬¡¢ËÙ¶È¡¢Î»ÖÃ¡¢°²×°Îó²î½Ç¡¢ÍÓÂÝÒÇÁãÆ«¡¢¼ÓËÙ¶È¼ÆÁãÆ«
        N = 18;
        H = zeros(9,N);
        H(1:3,4:6) = eye(3);
        H(4:6,7:9) = eye(3);
        H(7:9,1:3) = eye(3);
        P = diag([[1,1,1]*P0_phi, [1,1,1]*P0_dv, [P0_dlat,P0_dlon,P0_dh], [1,1,1]*P0_xi, [1,1,1]*P0_e, [1,1,1]*P0_d])^2;
        Q = diag([[1,1,1]*sigma_gyro*0.707, [1,1,1]*sigma_acc*0.707, [1/6378137,1/6378137,1]*sigma_acc*0.707*dt/2, [1,1,1]*0, [1,1,1]*0, [1,1,1]*0])^2 * dt^2;
        R = diag([[1,1,1]*sigma_v, [sigma_lat,sigma_lon,sigma_h], [1,1,1]*0.001])^2;
end
P0 = sqrt(diag(P))';

%% 3.Transfer alignment
n = (size(imu,1)-1)/2; %the number of inertial solving
nav = zeros(n,10); %[t, lat, lon, alt, vn, ve, vd, yaw, pitch, roll]
assembly_esti = zeros(n,4);
bias_esti = zeros(n,7);
filter_P = zeros(n,N);

p = traj(1,2:4) + p0_error; %deg
v = traj(1,5:7) + v0_error; %m/s
att = traj(1,8:10) + att0_error; %deg
pva0 = [p, v, att]; %record initial value
p(1:2) = p(1:2)/180*pi; %rad
att = att/180*pi; %rad
q = angle2quat(att(1), att(2), att(3));
qvp = [q, v, p]'; %column vector

dgyro = [0; 0; 0]; %gyro compensation, rad/s
dacc = [0; 0; 0]; %accelerometer compensation, m/s^2
assembly = [0, 0, 0]; %rad

X = zeros(N,1);

for k=1:n
    kj = 2*k+1;
    gyro0 = imu(kj-2, 2:4)'-dgyro; %rad/s
    gyro1 = imu(kj-1, 2:4)'-dgyro;
    gyro2 = imu(kj  , 2:4)'-dgyro;
    acc0  = imu(kj-2, 5:7)'-dacc; %m/s^2
    acc1  = imu(kj-1, 5:7)'-dacc;
    acc2  = imu(kj  , 5:7)'-dacc;
    
%     qvp = RK4(@ins_qvp_n, qvp, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
%     qvp(1:4) = quatnormalize(qvp(1:4)')';
    qvp = ins_qvp_n3(qvp, dt, [gyro0;acc0],[gyro1;acc1],[gyro2;acc2]);
    %---------------------------------------------------------------------%
    
    %--filtering
    Cab = angle2dcm(assembly(1), assembly(2), assembly(3));
    switch system_model
        case 'model_1'
            Phi = state_fun_1(qvp, acc1, dt);
            Z = qvp(5:10) - [measure(k+1,5:7)'; measure(k+1,2:4)'.*[pi/180;pi/180;1]];
        case 'model_2'
            Phi = state_fun_2(qvp, acc1, dt);
            Cna = angle2dcm(measure(k+1,8)/180*pi, measure(k+1,9)/180*pi, measure(k+1,10)/180*pi);
            Cbn = quat2dcm(qvp(1:4)')';
            Cna = Cab*Cna;
            C = Cbn*Cna;
            Z = [qvp(5:7)-measure(k+1,5:7)'; qvp(8:10)-measure(k+1,2:4)'.*[pi/180;pi/180;1]; C(2,3);C(3,1);C(1,2)];
            H(7:9,10:12) = measure_H(Cbn, Cna);
        case 'model_3'
            Phi = state_fun_3(qvp, acc1, dt);
            Cna = angle2dcm(measure(k+1,8)/180*pi, measure(k+1,9)/180*pi, measure(k+1,10)/180*pi);
            Cbn = quat2dcm(qvp(1:4)')';
            Cna = Cab*Cna;
            C = Cbn*Cna;
            Z = [qvp(5:7)-measure(k+1,5:7)'; qvp(8:10)-measure(k+1,2:4)'.*[pi/180;pi/180;1]; C(2,3);C(3,1);C(1,2)];
            H(7:9,10:12) = measure_H(Cbn, Cna);
        case 'model_4'
            Phi = state_fun_4(qvp, acc1, dt);
            Z = qvp(5:10) - [measure(k+1,5:7)'; measure(k+1,2:4)'.*[pi/180;pi/180;1]];
        case 'model_5'
            Phi = state_fun_5(qvp, acc1, dt);
            Cna = angle2dcm(measure(k+1,8)/180*pi, measure(k+1,9)/180*pi, measure(k+1,10)/180*pi);
            Cbn = quat2dcm(qvp(1:4)')';
            Cna = Cab*Cna;
            C = Cbn*Cna;
            Z = [qvp(5:7)-measure(k+1,5:7)'; qvp(8:10)-measure(k+1,2:4)'.*[pi/180;pi/180;1]; C(2,3);C(3,1);C(1,2)];
            H(7:9,10:12) = measure_H(Cbn, Cna);
    end
    [X, P] = kalman_filter(Phi, X, P, Q, H, Z, R);
    filter_P(k,:) = sqrt(diag(P))';
    
    %--adjust
    if norm(X(1:3))>0
        phi = norm(X(1:3));
        qc = [cos(phi/2), X(1:3)'/phi*sin(phi/2)];
        qvp(1:4) = quatmultiply(qc, qvp(1:4)')';
    end
    qvp(5:10) = qvp(5:10) - X(4:9);
    X(1:9) = zeros(9,1);
    switch system_model
        case 'model_2'
            if norm(X(10:12))>0
                phi = norm(X(10:12));
                qc = [cos(phi/2), X(10:12)'/phi*sin(phi/2)];
                Cab = quat2dcm(qc)*Cab;
                X(10:12) = zeros(3,1);
            end
        case 'model_3'
            if norm(X(10:12))>0
                phi = norm(X(10:12));
                qc = [cos(phi/2), X(10:12)'/phi*sin(phi/2)];
                Cab = quat2dcm(qc)*Cab;
                X(10:12) = zeros(3,1);
            end
            dgyro = dgyro + X(13:15);
            X(13:15) = zeros(3,1);
        case 'model_4'
            dgyro = dgyro + X(10:12);
            X(10:12) = zeros(3,1);
        case 'model_5'
            if norm(X(10:12))>0
                phi = norm(X(10:12));
                qc = [cos(phi/2), X(10:12)'/phi*sin(phi/2)];
                Cab = quat2dcm(qc)*Cab;
                X(10:12) = zeros(3,1);
            end
            dgyro = dgyro + X(13:15);
            dacc = dacc + X(16:18);
            X(13:18) = zeros(6,1);
    end
    [r1,r2,r3] = dcm2angle(Cab);
    assembly = [r1,r2,r3];
    assembly_esti(k,1) = k*dt;
    assembly_esti(k,2:4) = [r1,r2,r3] /pi*180;
    bias_esti(k,1) = k*dt;
    bias_esti(k,2:4) = dgyro' /pi*180;
    bias_esti(k,5:7) = dacc';
    
    %---------------------------------------------------------------------%
    nav(k,1) = k*dt;
    nav(k,2:3) = qvp(8:9)' /pi*180; %deg
    nav(k,4) = qvp(10); %m
    nav(k,5:7) = qvp(5:7)'; %m/s
    [r1,r2,r3] = quat2angle(qvp(1:4)');
    nav(k,8:10) = [r1,r2,r3] /pi*180; %deg
end
nav = [[0,pva0]; nav];
assembly_esti = [zeros(1,4); assembly_esti];
bias_esti = [zeros(1,7); bias_esti];
filter_P = [P0; filter_P];

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
% [filter_P(end,3)*3/pi*180*60; filter_P(end,2)*3/pi*180*60; filter_P(end,1)*3/pi*180*60;...
%  filter_P(end,15)*3/pi*180*3600; filter_P(end,16:18)'*3*100]

%% Function
function Phi = state_fun_1(qvp, fb, dt)
    global a f w
    v = qvp(5:7);
    lat = qvp(8);
    h = qvp(10);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    Cbn = quat2dcm(qvp(1:4)')';
    wien = [w*cos(lat); 0; -w*sin(lat)];
    wenn = [v(2)/Rn; -v(1)/Rm; -v(2)/Rn*tan(lat)];
    A = zeros(9);
    A(1:3,1:3) = -antisym(wien+wenn);
    A(4:6,1:3) = antisym(Cbn*fb);
    A(4:6,4:6) = -antisym(2*wien+wenn);
    A(7:9,4:6) = diag([1/Rm, sec(lat)/Rn, -1]);
    Phi = eye(9)+A*dt;%+(A*dt)^2/2;
end

function Phi = state_fun_2(qvp, fb, dt)
    global a f w
    v = qvp(5:7);
    lat = qvp(8);
    h = qvp(10);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    Cbn = quat2dcm(qvp(1:4)')';
    wien = [w*cos(lat); 0; -w*sin(lat)];
    wenn = [v(2)/Rn; -v(1)/Rm; -v(2)/Rn*tan(lat)];
    A = zeros(12);
    A(1:3,1:3) = -antisym(wien+wenn);
    A(4:6,1:3) = antisym(Cbn*fb);
    A(4:6,4:6) = -antisym(2*wien+wenn);
    A(7:9,4:6) = diag([1/Rm, sec(lat)/Rn, -1]);
    Phi = eye(12)+A*dt;%+(A*dt)^2/2;
end

function Phi = state_fun_3(qvp, fb, dt)
    global a f w
    v = qvp(5:7);
    lat = qvp(8);
    h = qvp(10);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    Cbn = quat2dcm(qvp(1:4)')';
    wien = [w*cos(lat); 0; -w*sin(lat)];
    wenn = [v(2)/Rn; -v(1)/Rm; -v(2)/Rn*tan(lat)];
    A = zeros(15);
    A(1:3,1:3) = -antisym(wien+wenn);
	A(1:3,13:15) = -Cbn;
    A(4:6,1:3) = antisym(Cbn*fb);
    A(4:6,4:6) = -antisym(2*wien+wenn);
    A(7:9,4:6) = diag([1/Rm, sec(lat)/Rn, -1]);
    Phi = eye(15)+A*dt;%+(A*dt)^2/2;
end

function Phi = state_fun_4(qvp, fb, dt)
    global a f w
    v = qvp(5:7);
    lat = qvp(8);
    h = qvp(10);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    Cbn = quat2dcm(qvp(1:4)')';
    wien = [w*cos(lat); 0; -w*sin(lat)];
    wenn = [v(2)/Rn; -v(1)/Rm; -v(2)/Rn*tan(lat)];
    A = zeros(12);
    A(1:3,1:3) = -antisym(wien+wenn);
	A(1:3,10:12) = -Cbn;
    A(4:6,1:3) = antisym(Cbn*fb);
    A(4:6,4:6) = -antisym(2*wien+wenn);
    A(7:9,4:6) = diag([1/Rm, sec(lat)/Rn, -1]);
    Phi = eye(12)+A*dt;%+(A*dt)^2/2;
end

function Phi = state_fun_5(qvp, fb, dt)
    global a f w
    v = qvp(5:7);
    lat = qvp(8);
    h = qvp(10);
    Rm = (1-f)^2*a / (1-(2-f)*f*sin(lat)^2)^1.5 + h;
    Rn =         a / (1-(2-f)*f*sin(lat)^2)^0.5 + h;
    Cbn = quat2dcm(qvp(1:4)')';
    wien = [w*cos(lat); 0; -w*sin(lat)];
    wenn = [v(2)/Rn; -v(1)/Rm; -v(2)/Rn*tan(lat)];
    A = zeros(18);
    A(1:3,1:3) = -antisym(wien+wenn);
	A(1:3,13:15) = -Cbn;
    A(4:6,1:3) = antisym(Cbn*fb);
    A(4:6,4:6) = -antisym(2*wien+wenn);
    A(4:6,16:18) = Cbn;
    A(7:9,4:6) = diag([1/Rm, sec(lat)/Rn, -1]);
    Phi = eye(18)+A*dt;%+(A*dt)^2/2;
end

function H = measure_H(Cbn, Cna)
    H = eye(3);
    H(1,1) = Cbn(2,3)*Cna(2,3) - Cbn(2,2)*Cna(3,3);
    H(1,2) = Cbn(2,1)*Cna(3,3) - Cbn(2,3)*Cna(1,3);
    H(1,3) = Cbn(2,2)*Cna(1,3) - Cbn(2,1)*Cna(2,3);
    H(2,1) = Cbn(3,3)*Cna(2,1) - Cbn(3,2)*Cna(3,1);
    H(2,2) = Cbn(3,1)*Cna(3,1) - Cbn(3,3)*Cna(1,1);
    H(2,3) = Cbn(3,2)*Cna(1,1) - Cbn(3,1)*Cna(2,1);
    H(3,1) = Cbn(1,3)*Cna(2,2) - Cbn(1,2)*Cna(3,2);
    H(3,2) = Cbn(1,1)*Cna(3,2) - Cbn(1,3)*Cna(1,2);
    H(3,3) = Cbn(1,2)*Cna(1,2) - Cbn(1,1)*Cna(2,2);
end