lat = 45;
w = 15;
g = 9.8;
dt = 0.02;

N = 6;
A = zeros(N);
A(1:3,1:3) = -antisym([cosd(lat),0,-sind(lat)]) *w/3600/180*pi;
A(4:6,1:3) = antisym([0,0,-g]);
if N==9
    A(1:3,7:8) = -[1,0; 0,0; 0,1];
    A(6,9) = 1;
end
Phi = eye(N)+A*dt+(A*dt)^2/2;

H = zeros(3,N);
H(1:3,4:6) = eye(3);

sigma_gyro = 0.2; %deg/h
sigma_acc = 0.6e-3; %g
sigma_v = 1e-3; %m/s

if N==9
Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1,...
          [1,1,1]*(sigma_acc*g)*0.707*1,...
          [1,1]*(sigma_gyro/100/3600/180*pi), sigma_acc/100].^2) *dt^2;
elseif N==6
Q = diag([[1,1,1]*(sigma_gyro/3600/180*pi)*0.707*1,...
          [1,1,1]*(sigma_acc*g)*0.707*1].^2) *dt^2;
end
R = diag([1,1,1]*sigma_v^2);

A = Phi';
B = H';
[X,L,G] = dare(A,B,Q,R);

sd = sqrt(diag(X))*3;
sd(1:3) = sd(1:3)/pi*180;
if N==9
    sd(7:8) = sd(7:8)/pi*180*3600;
    sd(9) = sd(9)/g*1000;
end