clear all
close all
clc

D=[0 0 0 0;
   0 1 0 0;
   0 0 1 0;
   0 0 0 1]

G=[0 0 0 0;
   1 0 0 0;
   1 0 0 0
   1 0 0 0]

L= D-G

A=-L

e = eig(A)

% ejercicio_5 examen CCRM Módulo 1

x0 = [10; 10; -5; -5; 2; 2; -15; -15]; %Condiciones iniciales
%x0 = [1; 1; 2; 2; 5; 5; -4; -4]; %Condiciones iniciales
dt = 0.01; %Incremento de tiempo
tvector = dt:dt:15; %Vector de tiempo


figure(1) %Estados x1 (posición) / tiempo

[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,1),'b','LineWidth',2);

ylim([-10 10])
xlim([0 11])
grid on
hold on


[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,3),'r','LineWidth',2);

ylim([-10 10])
xlim([0 11])
grid on
hold on


[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,5),'m');

ylim([-10 10])
xlim([0 11])
grid on
hold on


[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,7),'k','LineWidth',2);
title('Posición');
xlabel('Tiempo');
ylabel('x_{11},x_{21},x_{31},x_{41}');
ylim([-10 10])
xlim([0 11])
grid on
hold on

legend('n1','n2','n3','n4')

%------------------------------------------------------------------------

figure(2) %Estados x2 (Velocidad) / tiempo

[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,2),'b','LineWidth',2);

ylim([-10 10])
xlim([0 11])
grid on
hold on


[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,4),'r','LineWidth',2);

ylim([-10 10])
xlim([0 11])
grid on
hold on


[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,6),'m');

ylim([-10 10])
xlim([0 11])
grid on
hold on


[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(t,x(:,8),'k','LineWidth',2);
title('Velocidad');
xlabel('Tiempo');
ylabel('x_{12},x_{22},x_{32},x_{42}');
ylim([-10 10])
xlim([0 11])
grid on
hold on

legend('n1','n2','n3','n4')


%-------------------------------------------------------------------------

figure(3) %Diagramas de fase

subplot(3,3,1)
[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(x(:,1),x(:,3),'b','LineWidth',2);

ylabel('x_{21}');
ylim([0.0000 0.0044])
xlim([0.0000 0.0044])
grid on
hold on


subplot(3,3,4)
[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(x(:,1),x(:,5),'b','LineWidth',2);

ylabel('x_{31}');
ylim([0.0000 0.0044])
xlim([0.0000 0.0044])
grid on
hold on


subplot(3,3,7)
[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(x(:,1),x(:,7),'b','LineWidth',2);

xlabel('x_{11}');
ylabel('x_{41}');
ylim([0.0000 0.0044])
xlim([0.0000 0.0044])
grid on
hold on

%-------------------------------------------------------------------------

subplot(3,3,5)
[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(x(:,3),x(:,5),'b','LineWidth',2);

%ylabel('x_{31}');
ylim([0.0000 0.0044])
xlim([0.0000 0.0044])
grid on
hold on


subplot(3,3,8)
[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(x(:,3),x(:,7),'b','LineWidth',2);

xlabel('x_{21}');
%ylabel('x_{41}');
ylim([0.0000 0.0044])
xlim([0.0000 0.0044])
grid on
hold on

%-------------------------------------------------------------------------


subplot(3,3,9)
[t,x] = ode45(@(t,x)ejer5(t,x), tvector, x0);
plot(x(:,5),x(:,7),'b','LineWidth',2);

xlabel('x_{31}');
%ylabel('x_{41}');
ylim([0.0000 0.0044])
xlim([0.0000 0.0044])
grid on
hold on