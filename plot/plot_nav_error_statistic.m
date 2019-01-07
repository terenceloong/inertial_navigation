x = nav_error_1600s;

sd = sqrt(var(x,0,1))*3;

disp('pos_error (km)')
disp(sd(1:3)'/1000)
disp('vel_error (m/s)')
disp(sd(4:6)')
disp('att_error ('')')
disp(sd(7:9)'*60)

% figure
% subplot(3,1,1)
% histogram(x(:,1))
% title('导航位置误差统计直方图');
% xlabel('北向位置误差(m)');
% subplot(3,1,2)
% histogram(x(:,2))
% xlabel('东向位置误差(m)');
% subplot(3,1,3)
% histogram(x(:,3))
% xlabel('高度误差(m)');
% 
% figure
% subplot(3,1,1)
% histogram(x(:,4))
% title('导航速度误差统计直方图');
% xlabel('北向速度误差(m/s)');
% subplot(3,1,2)
% histogram(x(:,5))
% xlabel('东向速度误差(m/s)');
% subplot(3,1,3)
% histogram(x(:,6))
% xlabel('地向速度误差(m/s)');
% 
% figure
% subplot(3,1,1)
% histogram(x(:,7))
% title('导航姿态误差统计直方图');
% xlabel('航向角误差(\circ)');
% subplot(3,1,2)
% histogram(x(:,8))
% xlabel('俯仰角误差(\circ)');
% subplot(3,1,3)
% histogram(x(:,9))
% xlabel('滚转角误差(\circ)');

clearvars x