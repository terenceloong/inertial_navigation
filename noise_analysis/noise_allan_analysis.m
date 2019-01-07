fs = 100;
% y = noise_generator(3600*6, fs, 0, 0, 1, 100);
y = noise_generator(3600*6, fs, 0.2, 0.04, 0.0015, 100);
allan_variance(y, fs, 'omega');

%--100Hz
% 0.1 ---- 0.6, 0.2, 0.0055
% 0.03 --- 0.2. 0.04, 0.0015

%œ»≈‰N,K£¨‘Ÿ≈‰B