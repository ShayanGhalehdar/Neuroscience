dt = 0.01; % Simulation time step
Duration = 20000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
Cm = 1; % Membrane capacitance in micro Farads
gNa = 120; % in Siemens, maximum conductivity of Na+ Channel
gK = 36; % in Siemens, maximum conductivity of K+ Channel
gl = 0.3; % in Siemens, conductivity of leak Channel
ENa = 55; % in mv, Na+ nernst potential
EK = -72; % in mv, K+ nernst potential
El = -49.4; % in mv, nernst potential for leak channel
vRest = -60; % in mv, resting potential
V = vRest * ones(1,T); % Vector of output voltage
I = 10 * ones(1,T); % in uA, external stimulus (external current)
% for example: I(1:10000) = 2; % an input current pulse
n = zeros(1,T);
m = zeros(1,T);
h = zeros(1,T);


% uncomment the followong piece for minimum time ( 1.14ms ) --------------
%to see at least one spike                                       Q1_4
%                                                           --------------
%{                                                          
for i=114:T
    I(i)=0;
end
%}
% uncomment the followong piece for linear current          --------------
% from 6.323 to 110                                              Q1_7
%                                                           --------------
%
I1=6.323;
I2=110;
for i=1:T
    I(i)=I1+(I2-I1)*(i-1)/T;
end
%}
% uncomment each of followong pieces for other current forms--------------
%                                                                Q1_9
%                                                           --------------
%{
for i=1:T
    I(i)=10*sawtooth(i/500,0.5); %triangle
end
%}
%{
f=20;%Hz
I=10*sin(2*pi*f*t/1000); %sinusoidal %t(s)
%}
%{ 
%pulse
for i=1:ceil(T/3)
    I(i)=0;
end
for i=ceil(2*T/3):T
    I(i)=0;
end
%}
%{
I(1:T) = 10*chirp(t ,0.001,t(end),0.1); %chirp
%}

v=linspace(-80,10);
u = vRest - v;
alpha_n = (.1 * u + 1)./(exp(1 + .1 * u) - 1) / 10;
beta_n = .125 * exp(u/80);
alpha_m = (u+25) ./ (exp(2.5+.1*u)-1)/10;
beta_m = 4*exp(u/18);
alpha_h = .07 * exp(u/20);
beta_h = 1 ./ (1+exp(3 + .1*u));

%-------------------------------------------
% Q1_1
%-------------------------------------------

T_n = 1 ./ (alpha_n + beta_n);
T_m = 1 ./ (alpha_m + beta_m);
T_h = 1 ./ (alpha_h + beta_h);

n_inf = alpha_n .* T_n;
m_inf = alpha_m .* T_m;
h_inf = alpha_h .* T_h;
%{
figure
hold on
title('Time Constants')
xlabel('voltage')
ylabel('value')
plot(v,T_n,'y','LineWidth',1.2);
plot(v,T_m,'r','LineWidth',1.2);
plot(v,T_h,'b','LineWidth',1.2);
legend('T_n','T_m','T_h')
hold off

figure
hold on
title('Steady-State Values')
xlabel('voltage')
ylabel('value')
plot(v,n_inf,'y','LineWidth',1.2);
plot(v,m_inf,'r','LineWidth',1.2);
plot(v,h_inf,'b','LineWidth',1.2);
legend('n_{inf}','m_{inf}','h_{inf}')
hold off
%}

%-------------------------------------------
% Q1_2
%-------------------------------------------

V(1) = vRest;
for i = 1:T-1
    u = vRest-V(i);
    alpha_n = (.1 * u + 1)./(exp(1 + .1 * u) - 1) / 10;
    beta_n = .125 * exp(u/80);
    alpha_m = (u+25) ./ (exp(2.5+.1*u)-1)/10;
    beta_m = 4*exp(u/18);
    alpha_h = .07 * exp(u/20);
    beta_h = 1 ./ (1+exp(3 + .1*u));
    
    T_n = 1 / (alpha_n + beta_n);
    T_m = 1 / (alpha_m + beta_m);
    T_h = 1 / (alpha_h + beta_h);

    n_inf = alpha_n * T_n;
    m_inf = alpha_m * T_m;
    h_inf = alpha_h * T_h;
    
    if i==1 
        n(i) = n_inf;
        m(i) = m_inf;
        h(i) = h_inf;
    end
    
    dv = (dt/(-Cm))*(gl*(V(i)-El)+gK*n(i)^4*(V(i)-EK)+gNa*m(i)^3*h(i)*(V(i)-ENa)-I(i));
    dn = (dt/T_n)*(n_inf-n(i));
    dm = (dt/T_m)*(m_inf-m(i));
    dh = (dt/T_h)*(h_inf-h(i));
    
    V(i+1) = V(i)+dv;
    n(i+1) = n(i)+dn;
    m(i+1) = m(i)+dm;
    h(i+1) = h(i)+dh;

end
%{
figure
hold on
P1=subplot(2,1,1);
plot(t,V,'LineWidth',1.2);
title('Neuron Membrane Voltage');
xlabel('Time(ms)');
ylabel('Voltage(mV)');

P2=subplot(2,1,2);
plot(t,I,'LineWidth',1.2);
title('Neuron Membrane Current');
xlabel('Time(ms)');
ylabel('Current(uA)');
hold off
%}

%-------------------------------------------
% Q1_8
%-------------------------------------------
%{
figure
hold on
P1=subplot(3,1,1);
plot(n,V,'LineWidth',1.2);
title('n - Membrane Voltage Graph');
xlabel(P1, 'Voltage(mV)');
ylabel(P1, 'n');

P2=subplot(3,1,2);
plot(m,V,'LineWidth',1.2);
title('m - Membrane Voltage Graph');
xlabel(P2, 'Voltage(mV)');
ylabel(P2, 'm');

P2=subplot(3,1,3);
plot(h,V,'LineWidth',1.2);
title('h - Membrane Voltage Graph');
xlabel(P2, 'Voltage(mV)');
ylabel(P2, 'h');
hold off
%}

%-------------------------------------------
% Q1_10
%-------------------------------------------
F = zeros(1,T);
spikeActive=0;
s=0;
for i=1:T
    if V(i)>-35 && spikeActive==0
        s=s+1;
        spikeActive=1;
    end
    if V(i)<-35
        spikeActive=0;
    end
    F(i)=s;
end

figure
hold on
plot(I,F,'LineWidth',1.2);
title('F - I Graph');
xlabel('Current(uA)');
ylabel('Firing rate');
hold off