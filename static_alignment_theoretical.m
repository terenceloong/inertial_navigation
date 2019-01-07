X = fsolve(@fun, [0;0;0;0;0;0]);
% X = [psi;theta;gamma;en;ed;ad]

X(1:3) = X(1:3)/pi*180; %deg
X(6) = X(6)*1000; %mg
X

function y = fun(x)

lat = 45;
Cnb = angle2dcm(0/180*pi, 0/180*pi, 0/180*pi);
wien = [15*cosd(lat);0;-15*sind(lat)];
gn = [0;0;-1];
dgyro = [1; 0.1; -1];
dacc = [1; 1; 0.5]*1e-3;

y = zeros(6,1);
y(1:3) = angle2dcm(x(1),x(2),x(3))'*(Cnb*gn+dacc) - (gn+[0;0;x(6)]);
y(4:6) = angle2dcm(x(1),x(2),x(3))'*(Cnb*wien+dgyro) - (wien+[x(4);0;x(5)]);

end