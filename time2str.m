function str = time2str()

t = clock;
str = [num2str(t(1)),'-',num2str(t(2),'%02d'),'-',num2str(t(3),'%02d'),'_',...
       num2str(t(4),'%02d'),'-',num2str(t(5),'%02d'),'-',num2str(floor(t(6)),'%02d')];

end