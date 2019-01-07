function [y, b] = imu_error(x, fs, bias, N, B, K, type)

[n, m] = size(x);

% noise
noise = randn(n,m) *sqrt(fs)*N;

% instability
cn = dsp.ColoredNoise('Color','pink', 'SamplesPerFrame',n, 'NumChannels',m);
instability = cn()*B;

% rate random walk
% rate_rw = zeros(n,m);
% in = randn(n,m) /sqrt(fs)*K;
% for k=2:n
%     rate_rw(k,:) = rate_rw(k-1,:) + in(k,:);
% end
rate_rw = randn(n,m) /sqrt(fs)*K;
rate_rw = cumsum(rate_rw, 1);

% bias
b = ones(n,1)*bias + instability + rate_rw;

% output
if strcmp(type,'rate')
    y = x + b+noise;
else
    y = x + (b+noise)/fs;
end

end