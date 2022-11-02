close all;clc;clear all;
figure;hold on;

% This is just to simulate the system


%Inestable
A=[-1 2;2.2 1.7];
B=[2;1.6];

T = 15;

%K1 = [1.18188  -0.66629] %K PI online

K1 = [-0.35548  0.98831] %K PI online

x(:,1)=[0.5;-0.4];

for i=1:T
    x(:,i+1)=A*x(:,i)-B*K1*x(:,i);
end

K2 = [-0.36083  0.98709] %DLQR off line

x2(:,1)=[0.5;-0.4];

for i=1:T
    x2(:,i+1)=A*x2(:,i)-B*K2*x2(:,i);
end


figure(1)
hold on;
grid on;
plot( (1:T), x(2,1:T),"marker", "v", "markerEdgeColor", "k", ... 
     "markersize", 6, "linewidth", 2, "color","r");
     
plot( (1:T), x2(2,1:T),"marker", "s", "markerEdgeColor", "k", ... 
     "markersize", 6, "linewidth", 2, "linestyle", "--", "color","blue");

     
xlabel("tiempo [ k ]")
xlim ([0,T])
     
title({"Sistema inestable", "Muestras", T})
legend("PI - en línea: Estado x_{2}","DLQR: Estado x_{2}")
