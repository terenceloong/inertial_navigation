file_property = dir('imu.bin');
n = file_property.bytes/8/7;
imu = zeros(n,7);
fileID = fopen('imu.bin','r');
for k=1:n
    imu(k,:) = fread(fileID,7,'double')';
end
fclose(fileID);

file_property = dir('measure.bin');
n = file_property.bytes/8/10;
measure = zeros(n,10);
fileID = fopen('measure.bin','r');
for k=1:n
    measure(k,:) = fread(fileID,10,'double')';
end
fclose(fileID);

clearvars -except imu measure