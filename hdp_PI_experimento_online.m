close all;
clc;
clear all;

A=[0.9 0.1;0 0.8];
B=[0;1];

Q=diag ( [ 1 1 ] ) ;
R=2;
T = 700;   %Tiempo sobre el cual se ejecuta el algoritmo onlin

P=[1 0;0 1];
K=-inv(R+B'*P*B)*B'*P*A %Cálculo de K inicial
%K =[2, 0];

[ K0 , P0 ] = dlqr ( A , B , Q , R );   

N = 18; %Tamaño de la trayectoria
p = 3; %Número de terminos independientes

phi(1:p,1:N)=0; %phi(x_{k})
phi1(1:p,1:N)=0; %phi(x_{k+1})
Y(1:N,1)=0; %Recompensa

l=0;
con = 0;
ini = 1;
ini2 = 2;


[inp,wk]=michirp(0.1,3000,T+1,0.05);

figure(1)
plot(inp(:,2), "marker", "v", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","m");
title("Perturbación")
grid on;

K_it(1,fix(T/N))=0;

fix(T/N)

x(:,1)=[0.5; -1]; %condiciones iniciales


for i=1:T

u=K*x(:,i);
#u=K*x(:,i) + inp(i,2)*0.05;
x(:,i+1)=A*x(:,i)+B*u;
        
    Y(1,1)=Y(2,1);
    Y(2,1)=Y(3,1);
    Y(3,1)=Y(4,1);
    Y(4,1)=Y(5,1);
    Y(5,1)=Y(6,1);
    Y(6,1)=Y(7,1);
    Y(7,1)=Y(8,1);
    Y(8,1)=Y(9,1);
    Y(9,1)=Y(10,1);
    Y(10,1)=Y(11,1);
    Y(11,1)=Y(12,1);
    Y(12,1)=Y(13,1);
    Y(13,1)=Y(14,1);
    Y(14,1)=Y(15,1);
    Y(15,1)=Y(16,1);
    Y(16,1)=Y(17,1);
    Y(17,1)=Y(18,1);
    Y(18,1)=x(:,i)'*Q*x(:,i)+u'*R*u;
    
    phi(:,1)=phi(:,2);
    phi(:,2)=phi(:,3);
    phi(:,3)=phi(:,4);
    phi(:,4)=phi(:,5);
    phi(:,5)=phi(:,6);
    phi(:,6)=phi(:,7);
    phi(:,7)=phi(:,8);
    phi(:,8)=phi(:,9);
    phi(:,9)=phi(:,10);
    phi(:,10)=phi(:,11);
    phi(:,11)=phi(:,12);
    phi(:,12)=phi(:,13);
    phi(:,13)=phi(:,14);
    phi(:,14)=phi(:,15); 
    phi(:,15)=phi(:,16);
    phi(:,16)=phi(:,17);
    phi(:,17)=phi(:,18);
    phi(:,18)=[x(1,i)^2;2*x(1,i)*x(2,i);x(2,i)^2];
    
    phi1(:,1)=phi1(:,2);
    phi1(:,2)=phi1(:,3);
    phi1(:,3)=phi1(:,4);
    phi1(:,4)=phi1(:,5);
    phi1(:,5)=phi1(:,6);
    phi1(:,6)=phi1(:,7);
    phi1(:,7)=phi1(:,8);
    phi1(:,8)=phi1(:,9);
    phi1(:,9)=phi1(:,10);
    phi1(:,10)=phi1(:,11);
    phi1(:,11)=phi1(:,12);
    phi1(:,12)=phi1(:,13);
    phi1(:,13)=phi1(:,14);
    phi1(:,14)=phi1(:,15);
    phi1(:,15)=phi1(:,16);
    phi1(:,16)=phi1(:,17);
    phi1(:,17)=phi1(:,18);
    phi1(:,18)=[x(1,i+1)^2;2*x(1,i+1)*x(2,i+1);x(2,i+1)^2];
    
    if mod(i,N)==0
      
      con = con + 1;

      PHI = phi - phi1;
      %      DPHI = [ 2*x1(1,i) 0; 2*x1(2,i) 2*x1(1,i); 0 2*x1(2,i)];

      W=pinv(PHI*PHI')*PHI*Y;
      P=[W(1) W(2); W(2) W(3)];

      K_ = inv(R+B'*P*B)*B'*P*A; %Dinámica del sistema
      K = -K_;

      K_it(con) = norm(K0-K_);


      for i=1:2
        for j =1:2         
        P1(i, j+l) = P(i,j);         
        end
        P1;
      end
        
      l=l+2;
        
    end %if
    
end %for

K_
P

[ K0 , P0 ] = dlqr ( A , B , Q , R )


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
p21(:,1) = 0;
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



%Graficas

figure(2);
hold on;

plot( (1:T+1)-1, x(1,:),"marker", "v", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");
     
plot( (1:T+1)-1, x(2,:),"marker", "+", "markerEdgeColor", "k", ... 
      "markersize", 4, "linewidth", 2, "color","b");
xlabel("tiempo(k)")
%xlim ([0,43])
%ylim([-0.55, hold on;0.6])
title({"ADP (PI)"; T})
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
legend({"p11","p12", "p21", "p22"},  "location", "east")

line([1 16], [4.72993 4.72993], "linestyle", "--", "color", "r")
line([1 16], [0.67950 0.67950], "linestyle", "--", "color", "m")
line([1 16], [1.69146 1.69146], "linestyle", "--", "color", "black")
%xlim ([1,16])
grid on


figure(4);
hold on;

plot(K_it,"marker", "o", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");
xlabel("Iteración")
ylabel("|| K* - K||")
grid on