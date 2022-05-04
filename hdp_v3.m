%Este código es presentado en el segundo avance de tesis, como parte de la desmostración
%de la iteración de valor (HDP) usando identificación de parametros en linea sin PE.
%Como parte de la identificación se toma una estructura de aproximador paramétrico:
%Mínimos cuadrados por lotes. 

% Se toma como referencia del algoritmo crítico adaptativo:Landelius Dissertation and papers.


% Considerar 
% x(k+1)=[0      1; -0.16 -1]*x(k)+[0;1]*u(k)
% L=[.34 -2];
close all;
clear all;
clc;

A=[0.9 0.1;0 0.8];
B=[0;1];
%K=[0 1];
Q=diag ( [ 1 1 ] ) ;
R=2;
x=[0.5;-0.4]; %Condiciones iniciales
T = 200;
eig(A);


% Solo para simular el sistema
%P=1*eye(2,2)
P=[1 0;2 1];
K=-inv(R+B'*P*B)*B'*P*A

##%K=[0.2 1.8];
##%K=-K;
eig(A-B*K);

##for i=1:T
##    x(:,i+1)=A*x(:,i)+B*K*x(:,i);
##end

##x;
##figure(1); 
##hold on;
##plot( (1:T+1)-1, x(1,:),"marker", "v", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","r");
##     
##plot( (1:T+1)-1, x(2,:),"marker", "o", "markerEdgeColor", "k", ... 
##      "markersize", 4, "linewidth", 2, "color","b");
##%ylim([-0.55, 0.6])
##title({"Comportamiento del sistema con una ganancia K:"; K; "Iteraciones totales:";T})
##legend("Estado x_{1}","Estado x_{2}")
##grid on


N = 3; %Tamaño de trayectoria
p = 3; %Número de terminos independientes


phi(1:p,1:N)=0;
Y(1:N,1)=0;

l=0;
con = 0;
for i=1:T
    u=K*x(:,i);    
    x(:,i+1)=A*x(:,i)+B*u;
    
    Y(1,1)=Y(2,1);
    Y(2,1)=Y(3,1);
    Y(3,1)=x(:,i)'*Q*x(:,i)+u'*R*u+x(:,i+1)'*P*x(:,i+1);
    
    phi(:,1)=phi(:,2);
    phi(:,2)=phi(:,3);
    phi(:,3)=[x(1,i)^2;2*x(1,i)*x(2,i);x(2,i)^2];
    
    if mod(i,p)==0
        con = con + 1;
        C=phi*phi';
        rC(con) = rank(C); %Comprobando el rango de la matriz C
        q=phi*Y;        
        W=pinv(C)*q;
        P=[W(1) W(2);W(2) W(3)];
        K_ = inv(R+B'*P*B)*B'*P*A;
        K = -K_;
             
        for i=1:2
          for j =1:2         
          P1(i, j+l) = P(i,j);         
        end
        
      end
      l=l+2;
        
    end
end
P;
K_;
x;
l1 =1;
p11=zeros(1, con+1);
p11(:,1) = 1;
for j=2:con+1;
  p11(1, j)= P1(1,l1);
  l1=l1+2;
end

l2 =2;
p12=zeros(1, con+1);
p12(:,1) = 0;
for j=2:con+1;
  p12(1, j)= P1(1,l2);
  l2=l2+2;
end

l3 =1;
p21=zeros(1, con+1);
p21(:,1) = 2;
for j=2:con+1;
  p21(1, j)= P1(2,l3);
  l3=l3+2;
end


l4 =2;
p22=zeros(1, con+1);
p22(:,1) = 1;
for j=2:con+1;
  p22(1, j)= P1(2,l4);
  l4=l4+2;
end



[ K0 , P0 ] = dlqr ( A , B , Q , R ) 

figure(2);
hold on;

plot( (1:T+1)-1, x(1,:),"marker", "v", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");
     
plot( (1:T+1)-1, x(2,:),"marker", "o", "markerEdgeColor", "k", ... 
      "markersize", 4, "linewidth", 2, "color","b");
xlabel("tiempo(k)")
#xlim ([0,43])
ylim([-0.55, 0.6])
title({"HDP (VI)"; T})
legend("Estado x_{1}","Estado x_{2}")
grid on


figure(3);
hold on;
plot(p11,"marker", "o", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");

plot( p12,"marker", "o", "markerEdgeColor", "k", ... 
"markersize", 4, "linewidth", 2, "color","b");

plot( p21,"marker", "o", "markerEdgeColor", "k", ... 
"markersize", 4, "linewidth", 2, "color","m");

plot( p22,"marker", "o", "markerEdgeColor", "k", ... 
"markersize", 4, "linewidth", 2, "color","g");
xlabel("Iteración")
title("Convergencia de elementos de la matriz P")
legend({"p11","p12",  "p21", "p22"},  "location", "east")

line([1 16], [4.72993 4.72993], "linestyle", "--", "color", "r")
line([1 16], [0.67950 0.67950], "linestyle", "--", "color", "m")
line([1 16], [1.69146 1.69146], "linestyle", "--", "color", "black")
xlim ([1,16])
grid on


figure(4);
hold on;

plot(rC,"marker", "o", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");
ylim([1.5,3.5])
xlabel("Iteraciones")
title("Rango de la matriz phi(x_{k})*phi(x_{k})")
grid on

