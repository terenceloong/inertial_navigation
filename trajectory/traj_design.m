function cmd = traj_design(t)
cmd = zeros(1,6);

%%
% if 20<t && t<=30
%     cmd(2) = 1*(t-20);
% elseif 30<t && t<=40
%     cmd(2) = 10;
% elseif 40<t && t<=50
%     cmd(2) = 10-1*(t-40);
% end

% cmd(4) = 0;
% cmd(5) = 1;
% cmd(6) = 90;

%% 俯仰机动
% if 20<t
%     cmd(2) = t-20;
% end

%% 滚转机动
if 30<t && t<=36
    cmd(3) = 5*(t-30);
elseif 36<t && t<=44
    cmd(3) = 30;
elseif 44<t && t<=50
    cmd(3) = 30-5*(t-44);
end

%%
cmd(1:3) = cmd(1:3)/180*pi;
cmd(4:6) = [cmd(4)*cosd(cmd(6)), cmd(4)*sind(cmd(6)), cmd(5)];
end