clc
clear
close all
S = [2, 8, 3, 1];
L = length(S);
weights = cell(1, L-1);
weights{1} = [9, 3, -1;
    -2, 0.5, 1;
    -2, -0.5, 1;
    9, -3, -1;
    9, -3, 1;
    2, 0.5, 1;
    2, -0.5, 1;
    9, 3, 1];
weights{2} = [-3, 1, 0, 0, 1, 1, 0, 0, 1;
    1, 0, -1, -1, 0, 0, 0, 0, 0;
    1, 0, 0, 0, 0, 0, 1, 1, 0];
weights{3} = [-2, 1, 1, 1];
x = 0;
y = 2.1;
fprintf('X is %d and Y is %d \n',x,y);
Output = Perceptron(x, y, weights, S);
fprintf('Output is %d \n',Output);
rng=4;
delta=0.02;
X = -rng:delta:rng;
Y = -rng:delta:rng;
Z=zeros(length(X),length(Y));
for i=1:length(X)
    for j=1:length(Y)
        Z(j,i)= Perceptron(Y(i), X(j), weights, S);
    end
end
[x,y]=meshgrid(X,Y);
contourf(x,y, Z);
pcolor(x, y, Z);
shading interp
colorbar

function z = Perceptron(x, y, weights, S)
L = length(S);
m = size(x, 1); 
a = cell(1, L);
a{1} = [x'; y'];
for i = 2:L
    a{i-1} = [ones(1, m); a{i-1}];
    a{i} = sign(weights{i-1}*a{i-1});
end
z = a{L};
end