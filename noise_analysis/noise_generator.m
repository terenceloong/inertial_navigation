function y = noise_generator(T, fs, N, B, K, p)

n = T*fs;

% noise
noise = randn(n,1)*sqrt(fs)*N;

% instability
cn = dsp.ColoredNoise('Color','pink', 'SamplesPerFrame',n);
instability = cn()*B;

% rate random walk
% rate_rw  = zeros(n,1);
% in = randn(n,1)/sqrt(fs)*K;
% for k=2:n
%     rate_rw(k) = (1-1/fs/1000)*rate_rw(k-1) + in(k);
% end
rate_rw = randn(n,1)/sqrt(fs)*K;
rate_rw = cumsum(rate_rw);

% output
y = noise + instability + rate_rw;

% plot
if p~=0
    t = (0:p:(n-1))'/fs;
    figure
    plot(t, noise(1:p:end,1));
    hold on
    plot(t, instability(1:p:end,1));
    plot(t, rate_rw(1:p:end,1), 'LineWidth',3);
    set(gca, 'xlim', [t(1),t(end)])
    xlabel('\itt\rm(s)');
end

end