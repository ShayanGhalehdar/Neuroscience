%% Q2_1
clc
epsilon = 0.05;
trial=200;
r = zeros(1,trial);
v = zeros(1,trial);
for i=1:trial/2
    r(i)=1;
end
v(1)=0;
for i=1:trial-1
    v(i+1)=v(i)+epsilon*(r(i)-v(i));
end
tr=1:trial;
%{
figure
hold on
plot(tr,v);
xlabel("Trial")
ylabel("V")
hold off
%}

%% Q2_2
clc
epsilon = 0.05;
rewardProbability = 0.5;
trial=200;
r = zeros(1,trial);
v = zeros(1,trial);
for i=1:trial
    if rand()<rewardProbability
        r(i)=2;
    end
end
v(1)=0;
for i=1:trial-1
    v(i+1)=v(i)+epsilon*(r(i)-v(i));
end
tr=1:trial;
%{
figure
hold on
plot(tr,v);
xlabel("Trial")
ylabel("V")
hold off
%}

%% Q2_3
clc
epsilon = 0.05;
trial=200;
r = zeros(1,trial);
v1 = zeros(1,trial);
v2 = zeros(1,trial);
for i=1:trial
    r(i)=1;
end
v1(1)=1;
v2(1)=0;
for i=1:trial-1
    v1(i+1)=v1(i)+epsilon*(r(i)-v1(i)-v2(i));
    v2(i+1)=v2(i)+epsilon*(r(i)-v1(i)-v2(i));
end
tr=1:trial;
%{
figure
hold on
subplot(2,1,1)
plot(tr,v1);
xlabel("Trial")
ylabel("V1")

subplot(2,1,2)
plot(tr,v2);
xlabel("Trial")
ylabel("V2")
hold off
%}

%% Q2_4
clc
epsilon = 0.05;
trial=400;
r = zeros(1,trial);
v1 = zeros(1,trial);
v2 = zeros(1,trial);
v1(1)=0;
v2(1)=0;
for i=1:trial-1
    if rand()<0.5
        r=1;
        v1(i+1)=v1(i)+epsilon*(r-v1(i));
        v2(i+1)=v2(i);
    else
        r=0;
        v1(i+1)=v1(i)+epsilon*(r-v1(i)-v2(i));
        v2(i+1)=v2(i)+epsilon*(r-v1(i)-v2(i));
    end
end
tr=1:trial;
%
figure
hold on
subplot(2,1,1)
plot(tr,v1);
xlabel("Trial")
ylabel("V1")

subplot(2,1,2)
plot(tr,v2);
xlabel("Trial")
ylabel("V2")
hold off
%}

%% Q2_4
clc
epsilon = 0.05;
trial=200;
r = zeros(1,trial);
v1 = zeros(1,trial);
v2 = zeros(1,trial);
v1(1)=0;
v2(1)=0;
r=1;
for i=1:trial-1
    v1(i+1)=v1(i)+epsilon*(r-v1(i)-v2(i));
    v2(i+1)=v2(i)+epsilon*(r-v1(i)-v2(i));
end
tr=1:trial;
%
figure
hold on
subplot(2,1,1)
plot(tr,v1);
xlabel("Trial")
ylabel("V1")

subplot(2,1,2)
plot(tr,v2);
xlabel("Trial")
ylabel("V2")
hold off
%}