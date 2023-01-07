%% Q1 part 1-1&2
clc
clear
close all
x1=[1;1;1;1];
x2=[1;1;1;-1];
x3=[1;1;-1;-1];
W=x1*x1'+x2*x2'+x3*x3';
%part 2
s1=[1;-1;-1;1];
s2=sign_new(W*s1)
W*s2
s3=sign_new(W*s2)
s4=sign_new(W*s3)

%% Q1 part 2-1
clc
clear
close all
m = 6;
pictures = cell(1, m);
pictures{1} = 'train\astronaut.png';
pictures{2} = 'train\camera.png';
pictures{3} = 'train\chelsea.png';
pictures{4} = 'train\coffee.png';
pictures{5} = 'train\hourse.png';
pictures{6} = 'train\motorcycle.png';
train_pics = cell(1, m);
train_data = cell(1, m);
for i = 1:m
    train_pics{i} = imread(pictures{i});
    N = numel(train_pics{i});
    train_data{i} = (double(reshape(train_pics{i}, N, 1)))*2/255-1;
end
for i = 1:6
    subplot(2, 3, i);
    imagesc(train_pics{i});
    title('Train Data')        
end
% Q1 part 2-2
W = zeros(N,N);
for i = 1:m
    W = W + train_data{i}*(train_data{i})';
end
imagesc(W)
title('Weights')
colorbar()

% Q1 part 2-3
output = cell(1, m);
for i = 1:m
    output{i} = sign(W*train_data{i});
end
for i = 1:m
    subplot(2, 3, i)
    imagesc(reshape(output{i}, sqrt(N), sqrt(N)));
    title('Output')
end
error = zeros(1, m);
for i = 1:m
    error(i) = sum((train_data{i} ~= output{i}));
end
fprintf('Error Percent');
for i = 1:m
    fprintf('Error for %d: %0.2f pixels (%.2f%%)\n', i, error(i), error(i)/N*100);
end
% Q1 part 2-4
noisy_num = 3000;
noisy_data = cell(1, m);
for i = 1:m
    temp = train_data{i};
    index = randperm(N, noisy_num);
    temp(index) = temp(index)*(-1);
    noisy_data{i} = temp;
end
for i = 1:m
    subplot(2, 3, i)
    imagesc(reshape(noisy_data{i}, sqrt(N), sqrt(N)));
    title('Noisy Output')
end

% Q1 part 2-5
correlation = zeros(m, m);
for i = 1:m
    for j = 1:m
        correlation(i, j) = (train_data{i})'*(noisy_data{j})*(1/N);
    end
end
correlation
% Q1 part 2-6
output = cell(1, m);
for i = 1:m
    output{i} = sign(W*noisy_data{i});
end
for i = 1:m
    subplot(2, 3, i)
    imagesc(reshape(output{i}, sqrt(N), sqrt(N)));
    title('Output N=3000')
end
error = zeros(1, m);
for i = 1:m
    error(i) = sum((train_data{i} ~= output{i}));
end
fprintf('Error Percent\n');
for i = 1:m
    fprintf('Error for %d: %.2f pixels (%.2f%%)\n', i, error(i), error(i)/N)*100;
end
% Q1 part 2-7
noisy_num = 8000;
noisy_data = cell(1, m);
iterations = 20;
for i = 1:m
    temp = train_data{i};
    index = randperm(N, noisy_num);
    temp(index) = temp(index)*(-1);
    noisy_data{i} = temp;
end
output = noisy_data;
for iter = 1: iterations
    for i = 1:m
        output{i} = sign(W*output{i});
    end
end
for i = 1:m
    subplot(2, 3, i)
    imagesc(reshape(output{i}, sqrt(N), sqrt(N)));
    title('Output N=8000 iter=10')
end
error = zeros(1, m);
for i = 1:m
    error(i) = sum((train_data{i} ~= output{i}));
end
fprintf('Error Percent\n');
for i = 1:m
    fprintf('Error for %d: %.2f pixels (%.2f%%)\n', i, error(i), error(i)/N)*100;
end


function y=sign_new(x)
for i=1:length(x)
    if(x(i)>=0)
        y(i)=1;
    else
        y(i)=-1;
        
    end
end
y=y';
end