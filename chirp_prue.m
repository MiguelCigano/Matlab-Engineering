clear all;
close all;
clc;


T=200;
[inp,wk] = michirp(1,7000000,T+1,0.05);  
x = inp(:,2);
largo = length(x)
##
figure(1)
##plot( x,"marker", "*", "markerEdgeColor", "k", ... 
##      "markersize", 8, "linewidth", 1, "color","white");
plot(x,'--bs',...
    'LineWidth',1,...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor', 'magenta')
    
xlabel("Número de muestras")
hold on
xlim ([0, T])
grid on

figure(2)
hist(x,30, "facecolor", "yellow")
xlabel("Valores agrupados")
ylabel("Frecuencia")
grid on