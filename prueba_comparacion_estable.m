close all;clc;clear all;
figure;hold on;

% This is just to simulate the system

##%Estable
##A = [0.9 0.1;0 0.8];
##B = [0;1];

%Estable
A=[1.8980 -0.9084;1 0];
B=[1;0];

##%Estable (requiere PE)
##A=[1.1 -0.3; 1 0];
##B=[1; 0];

##%Estable (requiere PE)
##A=[0 1;-0.16 -1];
##B=[0;1];

##%Inestable
##A=[-1 2;2.2 1.7];
##B=[2;1.6];

T = 200;

%K1 = [1.18188  -0.66629] %K PI online

%K1 = [1.14529  -0.66850] %K PI online

%x(:,1)=[0.5;-0.4];

%for i=1:T
%    x(:,i+1)=A*x(:,i)-B*K1*x(:,i);
%end
[inp,wk]=michirp(0.01,690000,T+1,0.05);

K2 = [1.17051  -0.68524] %DLQR off line

x2(:,1)=[0.5;-0.4];

for i=1:T
    u=K2*x2(:,i) + inp(i,2)*0.102
    x2(:,i+1)=A*x2(:,i)-B*u;
end


%Q-learning

%K3 = [1.17051   -0.68523]
  
%x3(:,1)=[0.5;-0.4];

%for i=1:T
%    x3(:,i+1)=A*x3(:,i)-B*K3*x3(:,i);
%end
  
  
figure(1)
hold on;
grid on;

##plot( (1:T), x(1,1:T),"marker", "v", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","r");
    
plot( (1:T), x2(1,1:T),"marker", "s", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "linestyle", "-", "color","blue");

##plot( (1:T), x3(1,1:T),"marker", "d", "markerEdgeColor", "k", ... 
##     "markersize", 6, "linewidth", 2, "linestyle", ":", "color","black");

xlabel("Tiempo [ k ]")
xlim ([1,T])
     
title({"Sistema estable", "Muestras con PE", T})
%legend("PI - en línea: Estado x_{1}","DLQR: Estado x_{1}", "Q-learning: Estado x_{1}")
legend("Estado x_{1} - PE")


