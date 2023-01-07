%% Q2 part 1
clc
clear
close all
data = load('Q2_data.mat');
Fs = data.Fs;
trials = data.trials;
t = linspace(-50, 200, size(trials, 2));
figure()
raster_plot(trials(1, :), t)
title('Raster Plot only first')
figure()
raster_plot(trials(1:20, :), t)
title('Raster Plot 20 first trials')
figure()
raster_plot(trials, t)
title('Raster Plot all trials')
% Q2 part 3
window = 10;
PETH(trials, window, t);
title('PETH 10 ms')
%%
% Q2 part 4
figure()
window1 = 2.5;
PETH(trials, window1, t);
title('PETH window=5')
figure()
window2 = 10;
PETH(trials, window2, t);
title('PETH window=20')
figure()
window3 = 17.5;
PETH(trials, window3, t);
title('PETH window=35')


function raster_plot(data, t)
n = size(data, 1);
hold on
for i = 1:n
    idx = find(data(i, :) == 1);
    time = t(idx);
    spikes = zeros(size(time));
    plot([time; time], [spikes; spikes+1] + i-1, 'black');
end
plot([0;0], [0, n], 'LineStyle','--', 'LineWidth', 2,'color', 'b');
title('Raster Plot')
xlabel('t (s)')
ylabel('Trials')
end

function PETH(data, window, t)
n =size(data, 1);
t_init = t(1);
t_end = t(end);
edges = zeros(1, ceil((t_end-t_init)/window) + 1);
edges(1) = t_init - (window/2);
firing_rate = zeros(1,  ceil((t_end-t_init)/window));
count = 1;
while true
    edges(count+1) = min([t_end+(window/2), edges(count)+window]);
    index = find(t>=edges(count) & t<edges(count+1));
    firing_rate(count) = sum(sum(data(:, index)))*1000/window*(1/n);
    count = count + 1;
    if edges(count) == t_end+(window/2)
        break
    end
end
histogram('BinEdges',edges,'BinCounts',firing_rate, 'FaceAlpha', 1);
hold on
plot([0; 0], [0, 1.1*max(firing_rate)], 'LineStyle','--', 'LineWidth', 2,'color', 'r');
end