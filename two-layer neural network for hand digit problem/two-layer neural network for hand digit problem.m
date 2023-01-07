%% Q1 part 1
clear ; close all; clc
input_layer_size  = 400;  % 20x20 Input Images of Digits
hidden_layer_size = 25;   % 25 hidden units
num_labels = 10;          % 10 labels, from 1 to 10   
                          % (note that we have mapped "0" to label 10)
%  We start the exercise by first loading and visualizing the dataset. 
%  You will be working with a dataset that contains handwritten digits.

% Load Training Data
fprintf('Loading and Visualizing Data ...\n')

load('data.mat');
m = size(X, 1);

% Randomly select 100 data points to display
sel = randperm(size(X, 1));
sel = sel(1:100);

displayData(X(sel, :))

%% Q1 part2

data = load('data.mat');
x = data.X;
y = data.y;
m = size(x, 1);
n = size(x, 2);
preview = cell(1, 100);
idx = randperm(m, 100);
for i = 1:100
    preview{i} = reshape(x(idx(i), :), sqrt(n), sqrt(n));
end
x_train = zeros(300*10, n);
y_train = zeros(300*10, 1);

x_test = zeros(200*10, n);
y_test = zeros(200*10, 1);

for i = 1:10
    x_train((i-1)*300+1:i*300, :) = x((i-1)*500+1:(i-1)*500+300, :);
    y_train((i-1)*300+1:i*300) = y((i-1)*500+1:(i-1)*500+300);

    x_test((i-1)*200+1:i*200, :) = x((i-1)*500+301:i*500, :);
    y_test((i-1)*200+1:i*200) = y((i-1)*500+301:i*500);    
end


%% Q1 part3
S = [n, 25, 10];
lambda = 1;
L = length(S);
num_params = 0;
for i = 1:L-1
    num_params = num_params + (S(i)+1)*S(i+1);
end
weights = cell(1, L-1);
eps = 0.12;
for i = 1:L-1
    weights{i} = rand(S(i+1), S(i)+1)*2*eps - eps;
end
weights_vector = zeros(num_params, 1);
accum = 0;
for i = 1:L-1
    temp = weights{i};
    weights_vector(accum+1: accum+(S(i)+1)*S(i+1)) = temp(:);
    accum = accum +(S(i)+1)*S(i+1);
end

%% Q1 part 4,5,6
clc; close all;
input_layer_size  = 400;
hidden_layer_size = 25;   
num_labels = 10; 
costFunction = @(p) nnCostFunction(p, S, x_train, y_train, lambda, num_params,input_layer_size,hidden_layer_size,num_labels);

options = optimset('MaxIter', 200);
[weights_vector, cost] = fmincg(costFunction, weights_vector, options);


weights = cell(1, L-1);
accum = 0;
for i = 1:L-1
    weights{i} = reshape(weights_vector(accum+1:accum+S(i+1)*(S(i)+1)), S(i+1),S(i)+1);
    accum = accum + S(i+1)*(S(i)+1);
end



%% Q1 part 7
k = 5;
idx = randperm(size(x_test, 1), k*k);

pred_train = nnPredict(x_train, weights, L);
pred_test = nnPredict(x_test, weights, L);

figure()
title('Predictions of the model for test set')
for i = 1:k*k
    subplot(k, k, i);
    imshow(reshape(x_test(idx(i), :), sqrt(n), sqrt(n)));
    label = pred_test(idx(i));
    if label == 10
        label = 0;
    end
    title(num2str(label));
end
%% Q1 part 8
clc; close all;

second_layer = cell(1, S(2));
temp = weights{1};
for i = 1:S(2)
    second_layer{i} = reshape(temp(i, 2:end), sqrt(n), sqrt(n))+0.4;
end


montage(second_layer);
title('second layer')

acc_train = mean(pred_train == y_train);
acc_test = mean(pred_test == y_test);
fprintf('Accuracy of training set : %.3f%% \n', acc_train*100);
fprintf('Accuracy of test set :     %.3f%% \n', acc_test*100);


%% functions

function [h, display_array] = displayData(X, example_width)
if ~exist('example_width', 'var') || isempty(example_width) 
	example_width = round(sqrt(size(X, 2)));
end
colormap(gray);
[m n] = size(X);
example_height = (n / example_width);
display_rows = floor(sqrt(m));
display_cols = ceil(m / display_rows);
pad = 1;
display_array = - ones(pad + display_rows * (example_height + pad), ...
                       pad + display_cols * (example_width + pad));
curr_ex = 1;
for j = 1:display_rows
	for i = 1:display_cols
		if curr_ex > m, 
			break; 
        end
		max_val = max(abs(X(curr_ex, :)));
		display_array(pad + (j - 1) * (example_height + pad) + (1:example_height), ...
		              pad + (i - 1) * (example_width + pad) + (1:example_width)) = ...
						reshape(X(curr_ex, :), example_height, example_width) / max_val;
		curr_ex = curr_ex + 1;
	end
	if curr_ex > m, 
		break; 
	end
end
h = imagesc(display_array, [-1 1]);
axis image off
drawnow;
end

function [J, grad] = nnCostFunction(weights_vector , S, x, y, lambda, num_params,input_layer_size,hidden_layer_size,num_labels)
J = 0;
[m, ~] = size(x);
L = length(S);
a = cell(1, L);
delta = cell(1, L);
Delta = cell(1, L);
weights = cell(1, L-1);
accum = 0;
for i = 1:L-1
    weights{i} = reshape(weights_vector(accum+1:accum+S(i+1)*(S(i)+1)), S(i+1),S(i)+1);
    accum = accum + S(i+1)*(S(i)+1);
end
a{1} = x';
for i = 2:L
    a{i-1} = [ones(1, m);a{i-1}];
    a{i} = sigmoid(weights{i-1}*a{i-1});
end
labels = zeros(S(L), m);
labels((0:m-1)*S(L) + y') = 1;
delta{L} = a{L} - labels;
for i = L-1:-1:2
    temp = (weights{i})'*delta{i+1}.*a{i}.*(1 - a{i});
    delta{i} = temp(2:end, :); 
end
for i = 1:L-1
    Delta{i} = delta{i+1}*(a{i})'*(1/m);
    temp = weights{i};
    temp(:, 1) = 0;
    J = J + sum(sum(temp.^2))*lambda/(2*m);
    Delta{i} = Delta{i}+lambda/m*temp;
end
grad = zeros(num_params, 1);
accum = 0;
for i = 1:L-1
    temp = Delta{i};
    grad(accum+1: accum+(S(i)+1)*S(i+1)) = temp(:);
    accum = accum +(S(i)+1)*S(i+1);
end
J = J + (-1/m)*(sum(sum(labels.*log(a{L}) + (1-labels).*log(1-a{L}))));
end

function z = sigmoid(x)
    z = 1./(1+exp(-x));
end

function pred = nnPredict(x, weights, L)
m = size(x, 1);
a = cell(1, L);
a{1} = x';
for i = 2:L
    a{i-1} = [ones(1, m);a{i-1}];
    a{i} = sigmoid(weights{i-1}*a{i-1});
end
pred = a{L};
[~, pred] = max(pred);
pred = pred';
end