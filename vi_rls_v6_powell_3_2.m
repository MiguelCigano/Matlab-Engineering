close all;
clear all;
clc;

#{
Algoritmo de VI usando aproximadores en el crítico y en el actor para dinámica 
parcialmente conocida, (concimiento únicamente de la matriz B).
Asesores: Carlos Brizuelas, Alejandro Marquez.
Tesista: Jesus Miguel Martínez. 

Fecha: 5 Agosto del 2022

Código referente a iteración de valor online (HDP), usando Mínimos cuadrados recursivos (RLS).
Se considera el uso de Mínimos cuadrados ordinarios (LS) para la generación de theta de predicción (W_ant)
inicial a usar en el algoritmo de RLS.
Se modifica el codigo de vi_rls_v2_powell ya que no se estaba implementando 
correctamente el actor.
Para el actor se usará en este código decenso de gradiente.
#}


%Ecuaciones de estados de sistema----------------------------------------------

##%Estable
##A = [0.9 0.1;0 0.8];
##B = [0;1];

##%Estable
##A=[1.8980 -0.9084;1 0];
##B=[1;0];

##%Estable (requiere PE)
##A=[1.1 -0.3; 1 0];
##B=[1; 0];

##%Estable (requiere PE)
##A=[0 1;-0.16 -1];
##B=[0;1];

%Inestable
A=[-1 2;2.2 1.7];
B=[2;1.6];


%Criterio de minimización-------------------------------------------------------

Q = diag ( [ 1 1 ] );
R = 2;
eig(A)                                        %Polos de la matriz A en lazo abierto
gamma = 1;   
%Diseño de la ganancia inicial--------------------------------------------------

##[ K0 , P0 ] = dlqr( sqrt(gamma)*A , sqrt(gamma)*B , Q , R )
##
##P_e = P0*0.7
##K = -inv(R+B'*P_e*B)*B'*P_e*A
##K_1 = K(1, 1);
##K_2 = K(1, 2);

%Para el caso de la VI (iteración de valor) puede ser cualquier ganancia.


P_e = [1 0;0 1];
K = -inv(R+B'*P_e*B)*B'*P_e*A;
K_1 = K(1, 1);
K_2 = K(1, 2);

%Configuración de valores iniciales---------------------------------------------

T = 500;                                      %Tiempo total T
p = 3;                                        %Número de terminos independientes en nuestra matriz P_e
                                              %Números de experimentos
n_ls = 7;
n_ls2 = 3;                                     %Tamaño de lote para LS (mínimos cuadrados por lotes)
n_ls3 = 3;
xx = 1;
cont_ls_1 = 0;
l = 0;
con = 0;
exp = 0;
W = [P_e(1,1);P_e(1,2);P_e(2,2)];             %Pesos de W_j (estimación inicial) 
W_ant = [P_e(1,1);P_e(1,2);P_e(2,2)];         %Pesos de W_j+1 (predicción inicial)
##gamma = 1;                                    %Factor de descuento de la ecuación de Bellman
cont_2 = 0;
datos = 0;                                    %Número de iteraciones (muestras x_k) que se cumplen hasta la condición del error

phi= zeros(3,n_ls3);
Y = zeros(n_ls3,1);

%Datos del crítico--------------------------------------------------------------

f_olvido = 0.99;
a = 1 - f_olvido;
Iden = eye(p);
condicion_e_k = 10e-7;                       %Condición del error (error entre la estimación y la predicción)
convergencia = 1;                             %Convergencia = 0 , solo cuando se cumple la condición del error

%Datos del actor--------------------------------------------------------------
con_actor = 0;
beta=0.102;                                  %Parámetro de sintonización en el DG

%Señal de exitación-------------------------------------------------------------

[inp,wk]=michirp(0.1,690000,T+1,0.05);        %(f. inicial, f. final, muestras, amplitud)

%Condiciones iniciales----------------------------------------------------------

x(:,1)=[0.5;-0.4]; 

%Inicio del algoritmo-----------------------------------------------------------

   
  datos = 0;
  exp = exp+1;
  convergencia = 1;
  
  %P_n = 100*eye(3)
%Crítico------------------------------------------------------------------------

  for i=1:T
    
      %datos = datos + 1;
      
      u=K*x(:,i) + inp(i,2)*0.102;     
      %u=K*x(:,i);
      
      Y = circshift(Y,-1);
      phi = circshift(phi,-1, 2);
      
      x(:,i+1)= A*x(:,i) + B*u;
       
      r(i,:) = x(:,i)'*Q*x(:,i) + u'*R*u;
           
      phi(:,n_ls3) = [x(1,i)^2; x(1,i)*x(2,i) ; x(2,i)^2];
      %phi1(:,n_ls3)=[x(1,i+1)^2; 2*x(1,i+1)*x(2,i+1) ; x(2,i+1)^2];
      
      %Y(n_ls,:) = r(i,:) + gamma*W_ant'*phi1(:,n_ls);
      Y(n_ls3) = r(i,:) + gamma*x(:,i+1)'*P_e*x(:,i+1);
             
##      if i == n_ls
##        con_ls_1 = cont_ls_1 + 1;
##        
####        for k=1:n_ls
####          phi_(:,k) = phi(:,k);
####          Y_(k,:) = Y(k,:);
####        endfor 
##
##        C_ls = phi*phi';
##        W_ant = inv(C_ls)*phi*Y
##        theta_wls(:,n_ls) = W_ant;
##        P_n = phi*eye(n_ls)*phi';
##
##      end
    
                 
      if mod(i,n_ls3) == 0
     
        con = con + 1
        C_ls = phi*phi';
        rank(C_ls)
        W = inv(C_ls)*phi*Y;
        theta_wls(:,1) = W;
##        P_n = phi*eye(n_ls3)*phi';

      P_e=[theta_wls(1, 1) theta_wls(2,1)/2;theta_wls(2,1)/2 theta_wls(3,1)];
      K_e = inv(R+B'*P_e*B)*B'*P_e*A;
      K = -K_e;
        
##        %k=i;
##        
##        phi_k = phi;
##        r_k = r;       
##        phi1_k = phi1;
##        
##        for k=1:n_ls3
##        
##        phi_k = phi(:,k);
##        r_k = r(k,:);       
##        phi1_k = phi1(:,k);
##        
##        e_rls(:,k) = -W'*phi_k + r_k + gamma*W_ant'*phi1_k;
##        e_k = e_rls(:,k);
##        
####        if (abs(e_k)<condicion_e_k)
####          convergencia = 0;
####          datos
####        end  
##        
##        L1(:,k) = 1/f_olvido*P_n*phi_k*inv(1/a + phi_k'*1/f_olvido*P_n*phi_k);
##        L = L1(:, k);
##        
##        theta_wls(:, k) = W_ant + L*e_k;
##        P_n = 1/f_olvido*(Iden - L*phi_k')*P_n;
##        
##        W = theta_wls(:, k);
##        P_e = [theta_wls(1,k) theta_wls(2,k); theta_wls(2,k) theta_wls(3,k)];
##        
##        end
        %W_ant = W;  
      end
    
    #W_ant = W; 
    
  
  
##  if mod(i, n_ls3) == 0
##  i;
##  W;
##%Actor--------------------------------------------------------------------------
##        
##  cont_2 = cont_2 + 1;
####  U_j = K';
####  
####  for i=xx:n_ls2
####      con_actor = con_actor + 1;
####      sigma_k = x(:, i); %Sigma (x_k)
####      dphi(:,:) = [2*x(1,i+1) 0; 2*x(2,i+1) 2*x(1,i+1); 0 2*x(2,i+1)];
####
####      
####      diferencia = beta*sigma_k*(2*R*U_j'*sigma_k + gamma*B'*dphi'*W)';
####      U_j1  = U_j - diferencia;
####      K_es = U_j1;
####      U_j = K_es; %K_es ya está retroalimentada
####      
####      
####  endfor
####      xx = xx + n_ls;
####      n_ls2 = n_ls2 + n_ls;
####      
####      K_e = -K_es';
####      K = -K_e;
#### 
####  endif
##
##
      P_e=[theta_wls(1, n_ls3) theta_wls(2,n_ls3)/2;theta_wls(2,n_ls3)/2 theta_wls(3,n_ls3)];
      K_e = inv(R+B'*P_e*B)*B'*P_e*A;
      K = -K_e;
##
##  endif

##  if mod(i, n_ls3) == 0    
##      W_ant = W;
##  end 
      
  endfor
%-------------------------------------------------------------------------------
##      

##      
##        %Acomodar los valores de K para las gráficas
##        for c = 1:2
##          K1(1,c+l) = K_e(1,c);
##        end
##        
##        %Acomodar los valores de P para las gráficas
##        for f = 1:2
##          for c = 1:2         
##          P1(f,c+l) = P_e(f,c);
##          end       
##        end
##        
##        l=l+2;
    


K_e
P_e

[ K0 , P0 ] = dlqr( A , B , Q , R )

% Las siguientes lineas son para el armado de los vectores K_e y P_e------------
% K_e---------------------------------------------------------------------------
##
##g1 = 1;
##k11 = zeros(1, exp + 1); 
##k11(:, 1) = K_1;
##for j = 2:exp + 1
##  k11(1, j) = K1(1, g1);
##  g1 = g1 + 2;
##end
##
##g2 = 2;
##k12 = zeros(1, exp + 1);
##k12(:, 1) = K_2;
##for j = 2: exp + 1;
##  k12(1, j) = K1(1, g2);
##  g2 = g2 + 2;
##end
##
##
##
##% P_e---------------------------------------------------------------------------
##
##l1 = 1;
##p11 = zeros(1, exp + 1);
##p11(:,1) = 1; %P_11 ini
##for j=2:exp + 1;
##  p11(1, j)= P1(1,l1);
##  l1=l1+2;
##end
##
##l2 =2;
##p12=zeros(1, exp + 1);
##p12(:,1) = 0;
##for j=2: exp + 1;
##  p12(1, j)= P1(1,l2);
##  l2=l2+2;
##end
##
##l3 =1;
##p21=zeros(1, exp + 1);
##p21(:,1) = 0;
##for j=2:exp + 1;
##  p21(1, j)= P1(2,l3);
##  l3=l3+2;
##end
##
##l4 =2;
##p22=zeros(1, exp + 1);
##p22(:,1) = 1;
##for j=2: exp + 1;
##  p22(1, j)= P1(2,l4);
##  l4=l4+2;
##end
##
##%-------------------------------------------------------------------------------
##%Gráficas-----------------------------------------------------------------------
##
##figure(1);
##hold on;
##
##plot( (1:datos+1), x(1,1:datos+1),"marker", "v", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","r");
##     
##plot( (1:datos+1), x(2,1:datos+1),"marker", "o", "markerEdgeColor", "k", ... 
##      "markersize", 4, "linewidth", 2, "color","b");
##xlabel("tiempo(k)")
##xlim ([0, datos])
##
##title({"HDP (VI)"; T})
##legend("Estado x_{1}","Estado x_{2}")
##grid on
##
##figure(2);
##hold on;
##
##plot( (0:exp), p11,"marker", "o", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","r");
##
##plot( (0:exp), p12, "marker", "o", "markerEdgeColor", "k", ... 
##"markersize", 4, "linewidth", 2, "color","b");
##
##plot( (0:exp), p21, "marker", "o", "markerEdgeColor", "k", ... 
##"markersize", 4, "linewidth", 2, "color","m");
##
##plot( (0:exp), p22, "marker", "o", "markerEdgeColor", "k", ... 
##"markersize", 4, "linewidth", 2, "color","g");
##xlabel("Iteración")
##title("Convergencia de elementos de la matriz P")
##legend({"p11","p12",  "p21", "p22"},  "location", "east")
##
##line([0 exp], [ P0(1,1)  P0(1,1)], "linestyle", "--", "color", "r")
##line([0 exp], [P0(1,2) P0(1,2)], "linestyle", "--", "color", "m")
##line([0 exp], [P0(2,2) P0(2,2)], "linestyle", "--", "color", "black")
##xlim ([0,exp])
##grid on
##
##
##figure(3);
##hold on;
##
##plot( (0:(exp)), k11,"marker", "o", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","r");
##
##plot( (0:(exp)), k12, "marker", "o", "markerEdgeColor", "k", ... 
##"markersize", 4, "linewidth", 2, "color","b");
##
##xlabel("Iteración")
##title("Convergencia de elementos del vector de control K")
##legend({"k11","k12"},  "location", "east")
##
##line([0 num_exp], [ K0(1,1) K0(1,1) ], "linestyle", "--", "color", "r")
##line([0 num_exp], [K0(1,2) K0(1,2) ], "linestyle", "--", "color", "b")
##xlim ([0,exp])
##grid on
##
##
##
##
