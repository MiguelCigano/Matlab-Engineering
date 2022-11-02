
close all;
clear all;
clc;


 A=[1.8980 -0.9084;1 0];
 B=[1;0];

% A=[0.9 0.1;0 0.8];
% B=[0;1];

% %Estable (requiere PE)
% A=[0 1;-0.16 -1];
% B=[0;1];

##%Inestable
##A=[-1 2;2.2 1.7];
##B=[2;1.6];

Q=diag ( [ 1 1 ] );
R=2;
eig(A)

%P_e=[24 0.3;0.3 2];
P_e=[1 0;0 1];
K=-inv(R+B'*P_e*B)*B'*P_e*A


T = 250; %Tiempo total T
%N = 6; %Tamaño de trayectoria
p = 3; %Número de terminos independientes
n_ls = 7;
f_olvido = 0.99;
l = 0;
con = 0;

[inp,wk]=michirp(0.01,890000,T+1,0.05);

Gamma = eye(p,p)*1000; %varianza P

x(:,1)=[0.5;-0.4]; %Condiciones iniciales

for i=1:T
    u = K*x(:,i) + inp(i,2)*0.102;
    %u = K*x(:,i);
    x(:,i+1) = A*x(:,i)+B*u;
     
    
    Y(i,1) = x(:,i)'*Q*x(:,i)+u'*R*u+x(:,i+1)'*P_e*x(:,i+1);
    phi(:,i) = [x(1,i)^2;2*x(1,i)*x(2,i);x(2,i)^2];
    
    if i>=n_ls
     
      con = con + 1;
      if i == n_ls
        for k=1:n_ls
          phi_(:,k) = phi(:,k);
          Y_(k,:) = Y(k,:);
        end 

        C_minimos = phi_*phi_';
        W_ant = (pinv(C_minimos)*phi_)*Y_;
        W(:,n_ls) = W_ant;
      end
      
      W;
      j=i;
      L(:,j) = Gamma*phi(:,j)*inv(f_olvido + phi(:,j)'*Gamma*phi(:,j));
      W(:,j+1) = W(:,j) + L(:,j)*(Y(j) - phi(:,j)'*W(:,j));
      Gamma = 1/f_olvido * (Gamma - Gamma*inv(f_olvido + phi(:,j)'*Gamma*phi(:,j))*phi(:,j)*phi(:,j)'*Gamma);
      
    if mod(i,n_ls)==0
      P_e = [W(1, j+1) W(2,j+1);W(2,j+1) W(3,j+1)];
      K_e = inv(R+B'*P_e*B)*B'*P_e*A;
      K = -K_e;
    end
    
    
      for f = 1:2
        for c = 1:2         
        P1(f, c+l) = P_e(f,c);
        end       
      end
      l = l+2;
            
    end
    
    
end
K_e

P_e

[ K0 , P0 ] = dlqr( A , B , Q , R )

l1 =1;
p11=zeros(1, con+1);
p11(:,1) = 1; %P_11 ini
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

figure(3);
hold on;

plot( (0:T-(n_ls-1)), p11,"marker", "o", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");

plot( (0:T-(n_ls-1)), p12, "marker", "o", "markerEdgeColor", "k", ... 
"markersize", 4, "linewidth", 2, "color","b");

plot( (0:T-(n_ls-1)), p21, "marker", "o", "markerEdgeColor", "k", ... 
"markersize", 4, "linewidth", 2, "color","m");

plot( (0:T-(n_ls-1)), p22, "marker", "o", "markerEdgeColor", "k", ... 
"markersize", 4, "linewidth", 2, "color","g");
xlabel("Iteración")
title("Convergencia de elementos de la matriz P")
legend({"p11","p12",  "p21", "p22"},  "location", "east")

line([1 T], [ P0(1,1)  P0(1,1)], "linestyle", "--", "color", "r")
line([1 T], [P0(1,2) P0(1,2)], "linestyle", "--", "color", "m")
line([1 T], [P0(2,2) P0(2,2)], "linestyle", "--", "color", "black")
%xlim ([0,T+1])
grid on

figure(2);
hold on;

plot( (1:T+1)-1, x(1,:),"marker", "v", "markerEdgeColor", "k", ... 
     "markersize", 4, "linewidth", 2, "color","r");
     
plot( (1:T+1)-1, x(2,:),"marker", "o", "markerEdgeColor", "k", ... 
      "markersize", 4, "linewidth", 2, "color","b");
xlabel("tiempo(k)")
%xlim ([0,43])
%ylim([-0.55, 0.6])
title({"HDP (VI)"; T})
legend("Estado x_{1}","Estado x_{2}")
grid on