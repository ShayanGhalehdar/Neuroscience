%-------------------------------------------
% Q5_1
%-------------------------------------------
a=-10000:0.5:0;
x1=-sqrt(-a);
x2=sqrt(-a);
figure
hold on
title('Fixed Points - a')
xlabel('a')
ylabel('Fixed Points')
plot(a,x1,'r',a,x2,'r--','LineWidth',1.2);

