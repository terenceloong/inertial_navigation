N = 1000;

nav_error_200s = zeros(N,9);
nav_error_500s = zeros(N,9);
nav_error_800s = zeros(N,9);
nav_error_1200s = zeros(N,9);
nav_error_1600s = zeros(N,9);

for k=1:N
    clc;
    k
    
    imu = imu_delta;
    traj = traj_m;
    
    fs = 1/(imu(2,1)-imu(1,1));
    [imu(:,2:4), ~] = imu_error(imu(:,2:4), fs, [0,0,0], 0.006/3600/180*pi, 0.002/3600/180*pi, 0.000055/3600/180*pi, 'delta');
    [imu(:,5:7), ~] = imu_error(imu(:,5:7), fs, [0,0,0], 0.06*0.0098, 0.02*0.0098, 0.00055*0.0098, 'delta');
    %0.03--0.2, 0.04, 0.0015
    %0.1--0.6, 0.2, 0.0055
    
    p = traj(1,2:4);
    p(1:2) = p(1:2)/180*pi;
    v = traj(1,5:7);
    att = traj(1,8:10)/180*pi;
    q = angle2quat(att(1), att(2), att(3));
    avp = [q, v, p]';
    
    nav = test_inertial_solve(imu, avp, traj(:,4));
    
    index = find(nav(:,1)==200);
    nav_error_200s(k,:) = nav(index,2:end) - traj(index,2:end);
    index = find(nav(:,1)==500);
    nav_error_500s(k,:) = nav(index,2:end) - traj(index,2:end);
    index = find(nav(:,1)==800);
    nav_error_800s(k,:) = nav(index,2:end) - traj(index,2:end);
    index = find(nav(:,1)==1200);
    nav_error_1200s(k,:) = nav(index,2:end) - traj(index,2:end);
    index = find(nav(:,1)==1689);
    nav_error_1600s(k,:) = nav(index,2:end) - traj(index,2:end);
end

nav_error_200s = error_transform(nav_error_200s);
nav_error_500s = error_transform(nav_error_500s);
nav_error_800s = error_transform(nav_error_800s);
nav_error_1200s = error_transform(nav_error_1200s);
nav_error_1600s = error_transform(nav_error_1600s);

save(['./Data/statistic_inertial_solve_',time2str(),'.mat'],...
      'nav_error_200s','nav_error_500s','nav_error_800s','nav_error_1200s','nav_error_1600s');

function error = error_transform(error)
    error(:,1:2) = error(:,1:2)/180*pi*6378137; %lat and lon, deg to m
    for k=1:size(error,1)
        if error(k,7)>300
            error(k,7) = error(k,7)-360;
        elseif error(k,7)<-300
            error(k,7) = error(k,7)+360;
        end
        if error(k,9)>300
            error(k,9) = error(k,9)-360;
        elseif error(k,9)<-300
            error(k,9) = error(k,9)+360;
        end
    end
end