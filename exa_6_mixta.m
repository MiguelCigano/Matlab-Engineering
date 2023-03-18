clear all
close all
clc

G=[0 0 1 0 0 0 0;
   0 0 1 0 0 0 0;
   1 1 0 1 0 0 0; 
   1 0 1 0 0 0 1;
   1 0 0 0 0 0 0;
   1 0 0 0 0 0 0;
   1 0 0 1 0 0 0]

D=[1 0 0 0 0 0 0;
   0 1 0 0 0 0 0;
   0 0 3 0 0 0 0;
   0 0 0 3 0 0 0;
   0 0 0 0 1 0 0;
   0 0 0 0 0 1 0;
   0 0 0 0 0 0 2]

L= D-G

A=-L

e = eig(A)

% ejercicio_6 examen CCRM Módulo 1

x0 = [10; 10; -5; -5; 2; 2; -15; -15; 6; 6; 11; 11; -3; -3]; %Condiciones iniciales

dt = 0.01; %Incremento de tiempo
tvector = dt:dt:15; %Vector de tiempo


figure(1) %Estados x1 (posición) / tiempo

[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,1),'b','LineWidth',2);

ylim([-12 12])
xlim([0 9])
grid on
hold on

[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,3),'r','LineWidth',2);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,5),'m');

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,7),'c','linewidth',1.5);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,9),'g','linewidth',1.5);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,11),'y', 'linewidth',2);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,13),'k','LineWidth',2);
title('Posición');
xlabel('Tiempo');
ylabel('x_{11},x_{21},x_{31},x_{41}');
ylim([-12 12])
xlim([0 9])
grid on
hold on

legend('n1','n2','n3','n4','n5','n6','n7')


%------------------------------------------------------------------------

figure(2) %Estados x2 (Velocidad) / tiempo


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,2),'b','LineWidth',2);

ylim([-12 12])
xlim([0 9])
grid on
hold on

[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,4),'r','LineWidth',2);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,6),'m');

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,8),'c','linewidth',1.5);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,10),'g','linewidth',1.5);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,12),'y', 'linewidth',2);

ylim([-12 12])
xlim([0 9])
grid on
hold on


[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(t,x(:,14),'k','LineWidth',2);
title('Velocidad');
xlabel('Tiempo');
ylabel('x_{11},x_{21},x_{31},x_{41}');
ylim([-12 12])
xlim([0 9])
grid on
hold on

legend('n1','n2','n3','n4','n5','n6','n7')



%-------------------------------------------------------------------------

figure(3) %Diagramas de fase

subplot(3,6,1)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,1),x(:,3),'b','LineWidth',2);

ylabel('x_{21}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,7)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,1),x(:,5),'b','LineWidth',2);

ylabel('x_{31}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,13)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,1),x(:,7),'b','LineWidth',2);

ylabel('x_{41}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on

%-------------------------------------------------------------------------

figure(4)

subplot(3,6,1)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,1),x(:,9),'b','LineWidth',2);

ylabel('x_{51}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,7)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,1),x(:,11),'b','LineWidth',2);

ylabel('x_{61}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,13)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,1),x(:,13),'b','LineWidth',2);

xlabel('x_{11}');
ylabel('x_{71}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on

%-------------------------------------------------------------------------
figure (3)

subplot(3,6,8)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,3),x(:,5),'b','LineWidth',2);

ylabel('x_{31}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,14)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,3),x(:,7),'b','LineWidth',2);

ylabel('x_{41}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on

%-------------------------------------------------------------------------

figure (4)

subplot(3,6,2)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,3),x(:,9),'b','LineWidth',2);

ylabel('x_{51}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,8)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,3),x(:,11),'b','LineWidth',2);

ylabel('x_{61}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,14)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,3),x(:,13),'b','LineWidth',2);

ylabel('x_{71}');
xlabel('x_{21}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


%-------------------------------------------------------------------------

figure (3)

subplot(3,6,15)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,5),x(:,7),'b','LineWidth',2);

ylabel('x_{41}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on

%-------------------------------------------------------------------------

figure (4)

subplot(3,6,3)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,5),x(:,9),'b','LineWidth',2);

ylabel('x_{51}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,9)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,5),x(:,11),'b','LineWidth',2);

ylabel('x_{61}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,15)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,5),x(:,13),'b','LineWidth',2);

ylabel('x_{71}');
xlabel('x_{31}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on

%--------------

subplot(3,6,4)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,7),x(:,9),'b','LineWidth',2);

ylabel('x_{51}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,10)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,7),x(:,1),'b','LineWidth',2);

ylabel('x_{61}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,16)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,7),x(:,13),'b','LineWidth',2);

ylabel('x_{71}');
xlabel('x_{41}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,11)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,9),x(:,11),'b','LineWidth',2);

ylabel('x_{61}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,17)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,9),x(:,13),'b','LineWidth',2);

ylabel('x_{71}');
xlabel('x_{51}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on


subplot(3,6,18)
[t,x] = ode45(@(t,x)ejer6(t,x), tvector, x0);
plot(x(:,11),x(:,13),'b','LineWidth',2);

ylabel('x_{71}');
xlabel('x_{61}');
ylim([-0.00299 0.0000])
xlim([-0.00299 0.0000])
grid on
hold on