%-------------------------------------------
% Q2_3_1
%-------------------------------------------

dt = 0.1; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
a=0.02;
b=0.2;
c=-65;
d=2;
h=15;
vRest = -60;
u=b*vRest;
I = zeros(1,T);
V = zeros(1,T);
V(1)=vRest;
for i = 1:T-1
    if t(i)>10
        I(i)=h;
    end
    dv = dt*(0.04*V(i)^2 +5*V(i) +140 -u + I(i));
    V(i+1) = V(i)+dv;
    du = dt*a*(b*V(i)-u);
    u = u+du;
    
    if V(i+1)>30
        V(i+1)=c;
        u=u+d;
    end 
end

figure
hold on
P1=subplot(2,1,1);
plot(t,V,'LineWidth',1.2);
title('Tonic Spiking Voltage');
xlabel(P1, 'Time(ms)');
ylabel(P1, 'Voltage(mV)');

P2=subplot(2,1,2);
plot(t,I,'LineWidth',1.2);
title('Tonic Spiking Current');
xlabel(P2, 'Time(ms)');
ylabel(P2, 'Current(uA)');
hold off
%}

%-------------------------------------------
% Q2_3_2 & Q2_4
%-------------------------------------------
dt = 0.1; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
a=0.02;
b=0.25;
c=-65;
d=6;
h=1;
vRest = -60;
u=zeros(1,T);
u(1)=b*vRest;
I = zeros(1,T);
V = zeros(1,T);

V(1)=vRest;
for i = 1:T-1
    if t(i)>10
        I(i)=h;
    end
    dv = dt*(0.04*V(i)^2 +5*V(i) +140 -u(i) + I(i));
    V(i+1) = V(i)+dv;
    du = dt*a*(b*V(i)-u(i));
    u(i+1) = u(i)+du;
    
    if V(i+1)>30
        V(i+1)=c;
        u(i+1)=u(i+1)+d;
    end 
end

figure
hold on
P1=subplot(2,1,1);
plot(t,V,'LineWidth',1.2);
title('Phasic Spiking Voltage');
xlabel(P1, 'Time(ms)');
ylabel(P1, 'Voltage(mV)');

P2=subplot(2,1,2);
plot(t,I,'LineWidth',1.2);
title('Phasic Spiking Current');
xlabel(P2, 'Time(ms)');
ylabel(P2, 'Current(uA)');
hold off

figure
hold on
plot(u,V,'LineWidth',1.2);
title('Phasic Spiking Voltage based on parameter "u"');
xlabel('u(mV)');
ylabel('V(mV)');
hold off

%-------------------------------------------
% Q2_3_3
%-------------------------------------------
dt = 0.1; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
a=0.02;
b=0.2;
c=-50;
d=2;
h=15;
vRest = -60;
u=b*vRest;
I = zeros(1,T);
V = zeros(1,T);
V(1)=vRest;
for i = 1:T-1
    if t(i)>10
        I(i)=h;
    end
    dv = dt*(0.04*V(i)^2 +5*V(i) +140 -u + I(i));
    V(i+1) = V(i)+dv;
    du = dt*a*(b*V(i)-u);
    u = u+du;
    
    if V(i+1)>30
        V(i+1)=c;
        u=u+d;
    end 
end

figure
hold on
P1=subplot(2,1,1);
plot(t,V,'LineWidth',1.2);
title('Tonic Bursting Voltage');
xlabel(P1, 'Time(ms)');
ylabel(P1, 'Voltage(mV)');

P2=subplot(2,1,2);
plot(t,I,'LineWidth',1.2);
title('Tonic Bursting Current');
xlabel(P2, 'Time(ms)');
ylabel(P2, 'Current(uA)');
hold off
%}
%-------------------------------------------
% Q2_3_4
%-------------------------------------------
dt = 0.1; % Simulation time step
Duration = 100; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
a=0.02;
b=0.25;
c=-55;
d=0.05;
h=0.6;
vRest = -60;
u=0.27*vRest;
I = zeros(1,T);
V = zeros(1,T);
V(1)=vRest;
for i = 1:T-1
    if t(i)>10
        I(i)=h;
    end
    dv = dt*(0.04*V(i)^2 +5*V(i) +140 -u + I(i));
    V(i+1) = V(i)+dv;
    du = dt*a*(b*V(i)-u);
    u = u+du;
    
    if V(i+1)>30
        V(i+1)=c;
        u=u+d;
    end 
end

figure
hold on
P1=subplot(2,1,1);
plot(t,V,'LineWidth',1.2);
title('Phasic Bursting Voltage');
xlabel(P1, 'Time(ms)');
ylabel(P1, 'Voltage(mV)');

P2=subplot(2,1,2);
plot(t,I,'LineWidth',1.2);
title('Phasic Bursting Current');
xlabel(P2, 'Time(ms)');
ylabel(P2, 'Current(uA)');
hold off
%}

%-------------------------------------------
% Q2_3_5
%-------------------------------------------
dt = 0.1; % Simulation time step
Duration = 300; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
a=0.02;
b=0.2;
c=-55;
d=4;
h=10;
vRest = -60;
u=.4*vRest;
I = zeros(1,T);
V = zeros(1,T);
V(1)=vRest;
for i = 1:T-1
    if t(i)>10
        I(i)=h;
    end
    dv = dt*(0.04*V(i)^2 +5*V(i) +140 -u + I(i));
    V(i+1) = V(i)+dv;
    du = dt*a*(b*V(i)-u);
    u = u+du;
    
    if V(i+1)>30
        V(i+1)=c;
        u=u+d;
    end 
end

figure
hold on
P1=subplot(2,1,1);
plot(t,V,'LineWidth',1.2);
title('Mixed Model Voltage');
xlabel(P1, 'Time(ms)');
ylabel(P1, 'Voltage(mV)');

P2=subplot(2,1,2);
plot(t,I,'LineWidth',1.2);
title('Mixed Model Current');
xlabel(P2, 'Time(ms)');
ylabel(P2, 'Current(uA)');
hold off
%}