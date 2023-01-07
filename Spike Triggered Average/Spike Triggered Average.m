clc
clear variables
close all

load('Q3_data.mat')

t = 0:1/Fs:1;
figure
plot(t*1000, Stim(1:length(t)))
xlabel('t(ms)'); title('Stimulus Signal')

p = randperm(length(Spike_times), 20);
figure
for n = 1 : 20
    subplot(4, 5, n)
    t = Spike_times(p(n))-0.075:1/Fs:Spike_times(p(n));
    plot(t, Stim(int32(t*Fs)))
    title(strcat('Spike No. ', num2str(p(n))))
end

stim_pattern = zeros(length(Spike_times), 0.075*Fs+1);
for n = 1 : length(Spike_times)
    stim_pattern(n, :) = Stim(int32(((Spike_times(n)-0.075):1/Fs:Spike_times(n))*Fs));
end
spike_triggered_avg = mean(stim_pattern, 1);
figure
t = 0:1/Fs:0.075;
plot(t*1000, spike_triggered_avg)
title('Spike-Triggered Average')
xlabel('t(ms)')

figure
for n = 1 : 20
    plot(t*1000, Stim(int32((t+Spike_times(p(n))-0.075)*Fs)), 'LineWidth', 0.1, 'Color', 'g')
    hold on
end
plot(t*1000, spike_triggered_avg, 'LineWidth', 1.2, 'Color', 'r')
hold off
title('Spike-Triggered Average')
xlabel('t(ms)')
