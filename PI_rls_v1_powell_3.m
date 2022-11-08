close all;
clear all;
clc;

set(0, "defaultlinelinewidth", 2)

##05 falta imprimir las graficas de parametros
##04 de Nov Este código se deja listo para integrar las multigráficas

##Algoritmo de PI usando aproximadores en el crítico y en el actor para dinámica 
##parcialmente conocida, (concimiento únicamente de la matriz B).
##Asesores: Carlos Brizuelas, Alejandro Marquez.
##Tesista: Jesus Miguel Martínez. 
##
##Fecha: 29 Julio del 2022
##
##Código referente a iteración de política online (HDP), usando Mínimos cuadrados recursivos (RLS).
##Se considera el uso de Mínimos cuadrados ordinarios (LS) para la generación de theta de predicción (W)
##inicial a usar en el algoritmo de RLS.
##Para el actor se usará en este código decenso de gradiente.



%Ecuaciones de estado de sistema----------------------------------------------

##%Estable
##A = [0.9 0.1;0 0.8];
##B = [0;1];

%ESTABLE
A=[1.8980 -0.9084;1 0];
B=[1;0];

##%Estable (requiere PE)
##A=[1.1 -0.3;1 0];
##B=[1;0];

##%Estable (requiere PE)
##A=[0 1;-0.16 -1];
##B=[0;1];

%Inestable
##A=[-1 2;2.2 1.7];
##B=[2;1.6];




%Criterio de minimización-------------------------------------------------------

Q = diag ( [ 1 1 ] );
R = 2;
eig(A)                                        %Polos de la matriz A en lazo abierto

%Diseño de la ganancia inicial--------------------------------------------------


##[ K0 , P0 ] = dlqr( A , B , Q , R )

##P_e = P0*0.5
##K = -inv(R+B'*P_e*B)*B'*P_e*A
##K_1 = K(1, 1);
##K_2 = K(1, 2);

%P_e=[22 0.4;0.4 4];
P_e = [1 0;0 1]
K = -inv(R+B'*P_e*B)*B'*P_e*A
K_1 = -K(1, 1)
K_2 = -K(1, 2)

%Configuración de valores iniciales---------------------------------------------

T = 250;                                      %Tiempo total T
p = 3;                                        %Número de terminos independientes en nuestra matriz P_e
num_exp = 30;                                 %Números de experimentos
n_ls = 7;                                     %Tamaño de lote para LS (mínimos cuadrados por lotes)
l = 0;
con = 0;
exp = 0;
##W = [P_e(1,1);P_e(1,2);P_e(2,2)];             %Pesos de W_j (estimación inicial) 
##W_ant = [P_e(1,1);P_e(1,2);P_e(2,2)];         %Pesos de W_j+1 (predicción inicial)
gamma = 1;                                    %Factor de descuento de la ecuación de Bellman
cont_2 = 0;
datos = 0;                                    %Número de iteraciones (muestras x_k) que se cumplen hasta la condición del error


%Datos del crítico--------------------------------------------------------------

f_olvido = 0.99;
a = 1 - f_olvido;
Iden = eye(p);
condicion_e_k = 10e-6;                       %Condición del error (error entre la estimación y la predicción)
convergencia = 1;                             %Convergencia = 0 , solo cuando se cumple la condición del error
c_r = 1;

%Datos del actor--------------------------------------------------------------

beta=0.0102;                                  %Parámetro de sintonización en el DG

%Señal de exitación-------------------------------------------------------------

[inp,wk] = michirp(0.01,690000,T+1,0.05);        %(f. inicial, f. final, muestras, amplitud)

%Condiciones iniciales----------------------------------------------------------

x(:,1)=[0.5;-0.4]; 

%Inicio del algoritmo-----------------------------------------------------------

while exp<num_exp
  datos = 0;
  exp = exp+1
  convergencia = 1;
  dd = 1;

  P_ant = P_e;
  ip = 1;
  final_ruido = datos
%Crítico------------------------------------------------------------------------

  for i=1:T
    if (convergencia==1) %Convergencia a ayua a detener el algoritmo
      datos = datos + 1;
      final_ruido = datos;
      %u=K*x(:,i) + inp(i,2)*0.102;     
      %u=K*x(:,i);
      
      if exp >= 1
      u = K*x(:,i) + c_r*inp(i,2)*0.102;     
      end
      
##      if exp > 1
##        u=K*x(:,i);
##      end
      
      x(:,i+1) = A*x(:,i) + B*u;
       
      r(i,:) = x(:,i)'*Q*x(:,i) + u'*R*u;
           
      phi(:,i) = [x(1,i)^2; 2*x(1,i)*x(2,i) ; x(2,i)^2];
      phi1(:,i)=[x(1,i+1)^2; 2*x(1,i+1)*x(2,i+1) ; x(2,i+1)^2];
      
      Y(i,:) = r(i,:);
             
      if i == n_ls
        
        for k=1:n_ls
          phi_(:,k) = phi(:,k);
          phi1_k1_(:,k) = phi1(:,k);
          Y_(k,:) = Y(k,:);
        endfor 
        
        vec_act_ls_ = (phi_ - gamma*phi1_k1_);
        C_ls = vec_act_ls_*vec_act_ls_';
        W = inv(C_ls)*vec_act_ls_*Y_;
        theta_wls(:,n_ls) = W;
        P_n = vec_act_ls_*eye(n_ls)*vec_act_ls_';

      end
                 
      if i>n_ls
     
        con = con + 1;   
        
        k=i;
        
        phi_k = phi(:,k);
        r_k = r(k,:);       
        phi1_k1 = phi1(:,k);
        
        vec_act = (phi_k - gamma*phi1_k1);
        
        e_rls(:,k) = r_k - W'*(phi_k - gamma*phi1_k1); 
        e_k = e_rls(:,k);
        
        if (abs(e_k)<condicion_e_k)
          convergencia = 1;
          datos
          vec_datos(1,dd) = datos;
          c_r = 0;  %Detenemos la PE sobre la entrada de control
          dd = dd+1;
          final_ruido = min(vec_datos)
        end  
        
        L1(:,k) = 1/f_olvido*P_n*vec_act*inv(1/a + vec_act'*1/f_olvido*P_n*vec_act);
        L = L1(:,k);
        
        theta_wls(:, k) = W + L*e_k;
        P_n = 1/f_olvido*(Iden - L*vec_act')*P_n;
        
        W = theta_wls(:, k);
      
      end
    endif
  endfor
  
  
  x1 = x(1,:);
  x2 = x(2,:);  
  [nnn] = multigrafica (x1, x2, num_exp, exp, datos, final_ruido);

    
  P_e = [theta_wls(1,k) theta_wls(2,k); theta_wls(2,k) theta_wls(3,k)];
  
%Medición de la norma entre P_{j-1} y P_{j}-------------------------------------
  
  error_P(:,exp) = norm(P_e - P_ant);
  length(error_P);
  
%Actor--------------------------------------------------------------------------
        
  %cont_2 = cont_2 + 1;
  U_j = K';
  
  for i=1:datos
    
      sigma_k = x(:, i); %Sigma (x_k)
      dphi(:,:) = [2*x(1,i+1) 0; 2*x(2,i+1) 2*x(1,i+1); 0 2*x(2,i+1)];

      
      diferencia = beta*sigma_k*(2*R*U_j'*sigma_k + gamma*B'*dphi'*W)';
      U_j1  = U_j - diferencia;
      K_es = U_j1;
      U_j = K_es; %K_es ya está retroalimentada
      
  endfor
%-------------------------------------------------------------------------------
      error_K(:,exp) = norm(K' - K_es);
      K_e = -K_es';
      K = -K_e;

      
        %Acomodar los valores de K para las gráficas
        for c = 1:2
          K1(1,c+l) = K_e(1,c);
        end
        
        %Acomodar los valores de P para las gráficas
        for f = 1:2
          for c = 1:2         
          P1(f,c+l) = P_e(f,c);
          end       
        end
        
        l=l+2;
    
endwhile

##final_ruido = min(vec_datos)
e_rls;
K_e
P_e

[ K0 , P0 ] = dlqr( A , B , Q , R )

% Las siguientes lineas son para el armado de los vectores K_e y P_e------------
% K_e---------------------------------------------------------------------------

g1 = 1;
k11 = zeros(1, exp + 1); 
k11(:, 1) = K_1;
for j = 2:exp + 1
  k11(1, j) = K1(1, g1);
  g1 = g1 + 2;
end

g2 = 2;
k12 = zeros(1, exp + 1);
k12(:, 1) = K_2;
for j = 2: exp + 1;
  k12(1, j) = K1(1, g2);
  g2 = g2 + 2;
end



% P_e---------------------------------------------------------------------------

l1 = 1;
p11 = zeros(1, exp + 1);
p11(:,1) = 1; %P_11 ini
for j=2:exp + 1;
  p11(1, j)= P1(1,l1);
  l1=l1+2;
end

l2 =2;
p12=zeros(1, exp + 1);
p12(:,1) = 0;
for j=2: exp + 1;
  p12(1, j)= P1(1,l2);
  l2=l2+2;
end

l3 =1;
p21=zeros(1, exp + 1);
p21(:,1) = 0;
for j=2:exp + 1;
  p21(1, j)= P1(2,l3);
  l3=l3+2;
end

l4 =2;
p22=zeros(1, exp + 1);
p22(:,1) = 1;
for j=2: exp + 1;
  p22(1, j)= P1(2,l4);
  l4=l4+2;
end

%-------------------------------------------------------------------------------
%Gráficas-----------------------------------------------------------------------
##
##figure(1);
##hold on;
##
##
##plot( (1:datos+1), x(1,1:datos+1),"marker", "v", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","r");
##     
##plot( (1:datos+1), x(2,1:datos+1),"marker", "o", "markerEdgeColor", "k", ... 
##      "markersize", 4, "linewidth", 2, "color","b");
##xlabel("tiempo [ k ]")
##set(gca, 'FontSize', 17)
##%line([final_ruido final_ruido], [0 1], "linestyle", "--", "linewidth", 2,  "color", "black")
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
##xlabel("Experimento j")
##ylabel("Parámetros de P")
##set(gca, 'FontSize', 17)
##%title("Convergencia de elementos de la matriz P")
##legend({"p11","p12",  "p21", "p22"},  "location", "east")
##
##line([0 exp], [ P0(1,1)  P0(1,1)], "linestyle", "--", "color", "r")
##line([0 exp], [P0(1,2) P0(1,2)], "linestyle", "--", "color", "m")
##line([0 exp], [P0(2,2) P0(2,2)], "linestyle", "--", "color", "black")
##xlim ([0,exp])
##ylim ([-6,11])
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
##xlabel("Experimento j")
##ylabel("Parámetros de K")
##set(gca, 'FontSize', 17)
##%title("Convergencia de elementos del vector de control K")
##legend({"k11","k12"},  "location", "east")
##
##line([0 num_exp], [ K0(1,1) K0(1,1) ], "linestyle", "--", "color", "r")
##line([0 num_exp], [K0(1,2) K0(1,2) ], "linestyle", "--", "color", "b")
##xlim ([0,exp])
##grid on
##
##figure(4);
##hold on;
##
##plot( (1:exp), error_P,"marker", "o", "markerEdgeColor", "k", ... 
##     "markersize", 4, "linewidth", 2, "color","m");
##xlabel("Experimento j")
##ylabel("||P_{j+1} - P_{j}||")
##set(gca, 'FontSize', 17)
##%title("Cambio medido entre P_{j-1} y P_{j} durante el proceso de aprendizaje")
##%legend({"||P_{j-1} - P_{j}||"},  "location", "east")
##grid on
##
##figure(5);
##hold on;
##
##plot( (1:exp), error_K,"marker", "o", "markerEdgeColor", "red", ... 
##     "markersize", 4, "linewidth", 2, "color","cyan");
##xlabel("Experimento j")
##ylabel("||K_{j+1} - K_{j}||")
##ylim([0, 0.32])
##set(gca, 'FontSize', 17)
##%title("Cambio medido entre K_{j-1} y K_{j} durante el proceso de aprendizaje")
##%legend({"||K_{j-1} - K_{j}||"},  "location", "east")
##grid on
##
##