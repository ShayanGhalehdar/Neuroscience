clc;
%-------------------------------------------
% Q2_2
%-------------------------------------------
dt = 0.01; % Simulation time step
Duration = 1000; % Simulation length
T = ceil(Duration/dt);
t = (1:T) * dt; % Simulation time points in ms
El = -80;
Tm = 20;
IeRm = 25;
rm = 50;
Vth = -54;
Tp = 10;
K = t/Tp .* exp(1-t/Tp);
rho = zeros(1,T);
deltaPeriod=10000;
for i=deltaPeriod:deltaPeriod:T
    rho(i)=0.005;
end
g = conv(K,rho);
Ese = 0;
Esi = -80;
Es = Ese;
V = zeros(1,T);
V(1) = -80;
for i = 1:T-1
    if(V(i)<Vth)
        dv = (dt/Tm)*(-((V(i)-El)+g(i)*(V(i)-Es)*rm)+IeRm);
        V(i+1)=V(i)+dv;
    else
        dv = (dt/0.05)*(-V(i)+IeRm);
        V(i+1) = V(i)+dv;
        if V(i+1)>0
            V(i+1)=El;
        end
    end
end

%{
figure
hold on
P1=subplot(4,1,1);
plot(t,K,'LineWidth',1.2);
title('K');
P2=subplot(4,1,2);
plot(t,rho,'LineWidth',1.2);
title('rho');
P3=subplot(4,1,3);
plot(t,g(1:T),'LineWidth',1.2);
title('g');
P4=subplot(4,1,4);
plot(t,V,'LineWidth',1.2);
title('V');
hold off
%}

%-------------------------------------------
% Q2_3
%-------------------------------------------
El = -75;
V1 = zeros(1,T);
V2 = zeros(1,T);
V1(1) = -80;
V2(1) = -60;
rho1 = zeros(1,T);
rho2 = zeros(1,T);
g1 = conv(K,rho1);
g2 = conv(K,rho2);
excAmount = 0.0049;
%
for i=1:T-1
    if V1(i)<Vth && V2(i)<Vth
        dv1 = (dt/Tm)*(-((V1(i)-El)+g1(i)*(V1(i)-Es)*rm)+IeRm);
        V1(i+1)=V1(i)+dv1;
        dv2 = (dt/Tm)*(-((V2(i)-El)+g2(i)*(V2(i)-Es)*rm)+IeRm);
        V2(i+1)=V2(i)+dv2;
        
    elseif V1(i)<Vth && V2(i)>Vth
        dv1 = (dt/Tm)*(-((V1(i)-El)+g1(i)*(V1(i)-Es)*rm)+IeRm);
        V1(i+1)=V1(i)+dv1;
        dv2 = (dt/0.05)*(-V2(i)+IeRm);
        V2(i+1) = V2(i)+dv2;
        if V2(i+1)>0
            V2(i+1)=El;
            rho1(i)=excAmount;
            g1 = conv(K,rho1);
        end
        
    elseif V1(i)>Vth && V2(i)<Vth
        dv2 = (dt/Tm)*(-((V2(i)-El)+g2(i)*(V2(i)-Es)*rm)+IeRm);
        V2(i+1)=V2(i)+dv2;
        dv1 = (dt/0.05)*(-V1(i)+IeRm);
        V1(i+1) = V1(i)+dv1;
        
        if V1(i+1)>0
            V1(i+1)=El;
            rho2(i)=excAmount;
            g2 = conv(K,rho2);
        end
        
    else
        dv1 = (dt/0.05)*(-V1(i)+IeRm);
        V1(i+1) = V1(i)+dv1;
        dv2 = (dt/0.05)*(-V2(i)+IeRm);
        V2(i+1) = V2(i)+dv2;
        if V1(i+1)>0 || V2(i+1)>0
            V1(i+1)=El;
            rho2(i)=0.005;
            g2 = conv(K,rho2);
            V2(i+1)=El;
            rho1(i)=0.005;
            g1 = conv(K,rho1);
        end
        
    end
end
%}
figure
hold on
plot(t,V1,'r',t,V2,'b','LineWidth',1.2)
title('Inhibitory Synapses')
xlabel('time')
ylabel('Neurons Potential')
xlim([0 300])
legend('V1','V2')
hold off
