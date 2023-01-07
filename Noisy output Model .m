%-------------------------------------------
% Q3_1
%-------------------------------------------
beta=1;
D=500; %number of dots
range=5; gamma=0.5;%don't manipulate
x=linspace(-range,range,D);
F = beta*(1+tanh(gamma*x));
f=zeros(1,D);
for i=1:D-1
    f(i)=(F(i+1)-F(i))/(2*range/D);
end
%{
figure
hold on
P1=subplot(2,1,1);
plot(x,F,'LineWidth',1.2);
title('CDF');

P2=subplot(2,1,2);
plot(x,f,'LineWidth',1.2);
title('PDF');
hold off
%}

%-------------------------------------------
% Q3_2
%-------------------------------------------
dt = 0.01; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
Tau=1;
vTh = -45;
vRest = -70;
RI = 200;
V = zeros(1,T);
V(1)=vRest;
refractory=1;
th=zeros(1,ceil(Duration));
for i=1:ceil(Duration)
    while 1==1
        a=ceil(rand*D);
        if rand*f(D/2)<f(a)
            th(i)=vTh+x(a);
            break
        end
    end
end
j=1;
for i = 1:T-1
    if refractory==1
        V(i+1)=V(i)+0.1;
        if V(i+1)>th(j)
            j=j+1;
            refractory=0;
        end
    end
    if refractory==0
        dv = (dt/Tau)*(-V(i)+RI);
        V(i+1) = V(i)+dv;
        if V(i+1)>10
            V(i+1)=vRest;
            refractory=1;
        end
    end
end
%
figure
hold on
title('Neuron Membrane Voltage')
ylabel('voltage(mV)')
xlabel('Time(ms)')
plot(t,V,'LineWidth',1.2);
hold off
%}

%-------------------------------------------
% Q3_3
%-------------------------------------------
%{
figure
hold on
histogram(th(1:j-1))
title('Number of spikes in different voltage thresholds')
ylabel('Number of spikes')
xlabel('Voltage Threshold(mV)')
hold off
%}
%-------------------------------------------
% Q3_4
%-------------------------------------------

R = 1;
I = zeros(1,T);
for i=1:T
    I(i)=10+(i-1)*(10/(T-1));
end

V = zeros(1,T);
V(1)=vRest;
refractory=1;
j=1;
fr=zeros(1,T);
fires=0;
for i = 1:T-1
    if refractory==1
        V(i+1)=V(i)+1;
        if V(i+1)>th(j)
            if ceil(i/100000)==i/100000
                fires=0;
            end
            fires=fires+1;
            fr(i)=fr(i)+fires;
            j=j+1;
            refractory=0;
        end
    end
    if refractory==0
        dv = (dt/Tau)*(-V(i)+R*I(i));
        V(i+1) = V(i)+dv;
        if V(i+1)>10
            V(i+1)=vRest;
            refractory=1;
        end
    end
end

% you should increase the duration
%{
figure
hold on
title('F-I')
ylabel('F')
xlabel('I(uA)')
plot(I,fr,'LineWidth',1.2);
hold off
%}
