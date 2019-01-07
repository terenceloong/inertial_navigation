fileID = fopen('imu.bin','w+');
n = size(imu,1);
for k=1:n
    fwrite(fileID, imu(k,:), 'double');
end
fclose(fileID);

fileID = fopen('measure.bin','w+');
n = size(measure,1);
for k=1:n
    fwrite(fileID, measure(k,:), 'double');
end
fclose(fileID);